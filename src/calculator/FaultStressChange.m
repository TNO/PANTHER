classdef FaultStressChange
    % Fault stress change object. Calculates and stores 
    % shear and normal stress changes resulting from dP or dT. Controls the
    % calculation of stress for uniform vs depth-dependent P, T, or gamma
    % Loes Buijze 12 - 04 - 2023

    properties
        dsn double
        dtau double
    end

    properties (Dependent)
    end

    methods
        function self = FaultStressChange(size_y, size_t)
            % initialize StressChange class. Stress change contains changes
            % in shear and total normal stress
            % INPUT
            % size_y    length depth array
            % size_t    number of time steps
            self.dsn = zeros(size_y, size_t);
            self.dtau = zeros(size_y, size_t);
        end

        function self = calc_stress_changes(self, params, y, dx, pressure, temperature, load_case)
            % calculator for stress changes due to P and/or T.
            % controls the calculation for uniform and non-uniform pressure
            % and temperature profiles, non-uniform elastic parameters
            % INPUT
            % pressure, temperature     pressure, temperature objects
            % params                    input parameters for 1 ensemble member
            % y                         depth array w.r.t. y_mid
            % load_case                 P, T, or PT
            % check if pressure or temperature change is more than a single
            % value imposed in the reservoir. e.g. if there is diffusion
            [vary_P, vary_T] = self.variable_PT(pressure, temperature);
            % check if the dip is variable with depth
            [vary_dip] = self.variable_dip(params);
            if and(contains(load_case,'P'), vary_P) || and(contains(load_case,'T'), vary_T)
                vary_PT = 1;
            else
                vary_PT = 0;
            end
            GF = initialize_greens_functions(params, y, dx, vary_PT, vary_dip);
            % stress changes due to pressure changes
            if contains(load_case,'P')
              dsigma_P = self.calc_stress_changes_dP(params, y, pressure, GF, vary_PT, vary_dip );
            else 
              dsigma_P.dsn = zeros(size(self.dsn)); 
              dsigma_P.dtau = zeros(size(self.dtau)); 
            end    
            if contains(load_case,'T')
                dsigma_T = self.calc_stress_changes_dT(params, y, temperature, GF, vary_PT, vary_dip);
            else 
                dsigma_T.dsn = zeros(size(self.dsn)); 
                dsigma_T.dtau = zeros(size(self.dtau)); 
            end    
            % self.dsn = -self.dsn;
            self.dsn = -(dsigma_P.dsn + dsigma_T.dsn);
            self.dtau = dsigma_P.dtau + dsigma_T.dtau;
        end

        function dsigma = calc_stress_changes_dP(self, params, y, pressure, GF, vary_PT, vary_dip )
            % sets stress calculation for pressure changes
            gamma_P = params.get_gamma_P; 
            if and(~vary_dip, ~vary_PT)
                dsigma = self.get_stress_change_uniform(params, GF{1}, y, pressure, gamma_P,'P');
            else
                dsigma = self.get_stress_change_nonuniform( GF, y, size(pressure.dp_fault, 2), pressure, gamma_P, 'P');
            end
        end

        function dsigma = calc_stress_changes_dT(self, params, y, temperature, GF, vary_PT, vary_dip )
            % sets stress calculation for temperature changes
            gamma_T = params.get_gamma_T; 
            if and(~vary_dip, ~vary_PT)
                dsigma = self.get_stress_change_uniform(params, GF{1}, y, temperature, gamma_T, 'T');
            else
                dsigma = self.get_stress_change_nonuniform(GF, y, size(temperature.dT_fault, 2), temperature, gamma_T, 'T');
            end
        end

        function self = get_stress_change_uniform(self, params, GF, y, PT_change, gamma, load)
            % calculate the stress change for each timestep, for uniform P or T change in the reservoir blocks
            % find the pressure in the HW or FW compartment
            if strcmp(load, 'P') 
                for i = 1 : size(PT_change.dp_FW,2)
                    i_mid = floor((params.top_FW_i(y) + params.base_FW_i(y))/2);    
                    dP_FW = PT_change.dp_FW(i_mid, i);      % [MPa] take the pressure at the reservoir compartment center (it will be uniform)
                    i_mid = floor((params.top_HW_i(y) + params.base_HW_i(y))/2);    
                    dP_HW = PT_change.dp_HW(i_mid, i);      % [MPa] take the pressure at the reservoir compartment center (it will be uniform)
                    [self.dsn(:,i), self.dtau(:,i)] = self.get_stress_change_component(GF, dP_FW, dP_HW, gamma);   
                end
            else 
                for i = 1 : size(PT_change.dT_FW,2)
                    i_mid = floor((params.top_FW_i(y) + params.base_FW_i(y))/2);    
                    dT_FW = PT_change.dT_FW(i_mid, i);      % [MPa] take the pressure at the reservoir compartment center (it will be uniform)
                    i_mid = floor((params.top_HW_i(y) + params.base_HW_i(y))/2);    
                    dT_HW = PT_change.dT_HW(i_mid, i);      % [MPa] take the pressure at the reservoir compartment center (it will be uniform)
                    [self.dsn(:,i), self.dtau(:,i)] = self.get_stress_change_component(GF, dT_FW, dT_HW, gamma);   
                end
            end
        end

        function self = get_stress_change_nonuniform(self, GF, y, n_times, PT_change, gamma_PT, load_case)
            % calculate shear and normal stress, by adding contributions
            % from all depth increments, which may have varyiable P, T, or
            % gamma with depth 
            for i = 1 : n_times
                dsn_temp = zeros(length(y),1);      % array for adding contributions of gridded depth blocks withs varying P,T, or gamma
                dtau_temp = zeros(length(y),1);

                % refactor. element-wise multiplication in get stress
                % change component
                for j = 1 : length(y)
                    if strcmp(load_case,'P')
                        dPT_FW = PT_change.dp_FW(j,i);
                        dPT_HW = PT_change.dp_HW(j,i);
                    else 
                        dPT_FW = PT_change.dT_FW(j,i);
                        dPT_HW = PT_change.dT_HW(j,i);
                    end
                    if length(gamma_PT) == 1
                        gamma = gamma_PT;
                    elseif length(gamma_P) == length(y)
                        gamma = gamma_PT(j);
                    end
                    [dsn_j, dtau_j] = self.get_stress_change_component(GF{j}, dPT_FW, dPT_HW, gamma);
                    dsn_temp = dsn_temp + dsn_j;
                    dtau_temp = dtau_temp + dtau_j;
                end
                self.dsn(:,i) = dsn_temp;
                self.dtau(:,i) = dtau_temp;
            end
            % plot(y, PT_change.dp_HW(:,2), y, PT_change.dp_FW(:,2),y, dsn_temp)
        end

        function [dsn, dtau] = get_stress_change_component(~, GF, dPT_FW, dPT_HW, gamma)
            % compute the shear and normal stress change for separate
            % blocks
            dsn_FW = GF.Gnorm_FW * dPT_FW .* gamma / (2*pi);
            dtau_FW = GF.Gshear_FW * dPT_FW .* gamma / (2*pi);
            dsn_HW = GF.Gnorm_HW * dPT_HW .* gamma / (2*pi);
            dtau_HW = GF.Gshear_HW * dPT_HW .* gamma / (2*pi);
            dsn = dsn_FW + dsn_HW;
            dtau = dtau_FW + dtau_HW;
        end


        function plot_greens_f(~, greens_f, y)
            h1 = figure(1); clf(h1); hold on 
            gfs = properties(greens_f);
            for i = 1 : length(gfs)
                plot(y, greens_f.(gfs{i}), 'LineWidth',1.5);
            end
            legend(gfs);
        end

        % function to check if pressure change is uniform (consider moving
        % to pressure object)
        function check_uniform(~, pressure)
            unique_P = unique(pressure.dp_FW);
            if length(unique_P) > (size(pressure.dp_FW, 2) + 1)
                error('Pressure change in footwall is not uniform');
            end
            unique_P = unique(pressure.dp_HW);
            if length(unique_P) > (size(pressure.dp_HW, 2) + 1)
                error('Pressure change in hanging wall is not uniform');
            end
        end

        function [vary_P, vary_T] = variable_PT(~, pressure, temperature)
            % check if pressure or temperature are non uniform with y
            % return true if they vary . move to pressure object
            vary_P = 0;
            vary_T = 0;
            if length(unique(pressure.dp_FW)) > size(pressure.dp_FW,2) + 1
                vary_P = 1;
            elseif length(unique(pressure.dp_HW)) > size(pressure.dp_FW,2) + 1
                vary_P = 1;
            end
            if length(unique(temperature.dT_FW)) > size(temperature.dT_FW,2) + 1
                vary_T = 1;
            elseif length(unique(temperature.dT_HW)) > size(temperature.dT_FW,2) + 1
                vary_T = 1;
            end
        end

        function [vary_dip] = variable_dip(~, params)
            if or(length(params.dip) > 1, iscell(params.dip))
                vary_dip = 1;
            else
                vary_dip = 0;
            end
        end

        function [dscu] = get_scu_change(self, f_s, cohesion)
            % return Shear Capacity Utilization 
            % TODO: check if this works for depth-dependent f_s
            dscu = self.dtau ./ (self.dsn .* f_s + cohesion); 
        end
    end

end

            
%             greens_f_conv = GreensFunctions(y);              % initialize Green's functions
%             mid_y = y(251);
%             mid_x = mid_y/(tan(params.dip*pi/180));
%             greens_f_conv = greens_f_conv.green_FW(xeval, y, params.dip, slice_thick, slice_throw, 1e9, mid_x, mid_y);
%             for i = 1 : size(pressure.dp_FW,2)
%                 dsn_FW_conv(:,i) = conv(greens_f_test.Gnorm_FW, (pressure.dp_FW(:,i) ),'same')* params.get_gamma_P / (2*pi);
%             end
