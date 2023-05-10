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
        function self = FaultStressChange(pressure)
            % initialize StressChange class. Stress change contains changes
            % in shear and total normal stress
            self.dsn = zeros(size(pressure.dp_fault));
            self.dtau = zeros(size(pressure.dp_fault));
        end

        function self = calc_stress_changes(self, params, y, dx, pressure, load_case)
            % calculator for stress changes due to P and/or T.
            % controls the calculation for uniform and non-uniform pressure
            % and temperature profiles, non-uniform elastic parameters
            [vary_P, ~] = self.variable_PT(pressure, zeros(size(pressure)));
            if or(length(params.dip) > 1, iscell(params.dip))
                vary_dip = 1;
            else
                vary_dip = 0;
            end
            GF = self.initialize_greens_functions(params, y, dx, vary_P, vary_dip);
            if and(~vary_dip, ~vary_P)
                self = self.get_stress_change_uniform(params, GF{1}, y, pressure);
            else
                self = self.get_stress_change_nonuniform(params, GF, y, pressure);
            end
            self.dsn = -self.dsn;
        end

        function self = get_stress_change_uniform(self, params, GF, y, pressure)
            % calculate the stress change for each timestep, for uniform P or T change in the reservoir blocks
            for i = 1 : size(pressure.dp_FW,2)
                ip = find(pressure.dp_FW(:,i) ~= 0, 1);
                if ~isempty(ip)
                    dP_FW = pressure.dp_FW(ip, i);
                else
                    dP_FW = 0;
                end
                % i_mid = floor((params.top_HW_i(y) + params.base_HW_i(y))/2);    
                % dP_FW = pressure.dp_FW(i_mid, i);
                i_mid = floor((params.top_HW_i(y) + params.base_HW_i(y))/2);    
                dP_HW = pressure.dp_HW(i_mid, i);      % [MPa] take the pressure at the reservoir compartment center (it will be uniform)
                gamma = params.get_gamma_P; 
                [self.dsn(:,i), self.dtau(:,i)] = self.get_stress_change_component(GF, dP_FW, dP_HW, gamma);   
            end
        end

        function self = get_stress_change_nonuniform(self, params, GF, y, pressure)
            % calculate shear and normal stress, by adding contributions
            % from all depth increments, which may have varyiable P, T, or
            % gamma with depth
            gamma_P = params.get_gamma_P;
            for i = 1 : size(pressure.dp_fault, 2)
                dsn_temp = zeros(length(y),1);      % array for adding contributions of gridded depth blocks withs varying P,T, or gamma
                dtau_temp = zeros(length(y),1);
                for j = 1 : length(y)
                    dPT_FW = pressure.dp_FW(j,i);
                    dPT_HW = pressure.dp_HW(j,i);
                    if length(gamma_P) == 1
                        gamma = gamma_P;
                    elseif length(gamma_P) == length(y)
                        gamma = gamma_P(j);
                    end
                    [dsn_j, dtau_j] = self.get_stress_change_component(GF{j}, dPT_FW, dPT_HW, gamma);
                    dsn_temp = dsn_temp + dsn_j;
                    dtau_temp = dtau_temp + dtau_j;
                end
                self.dsn(:,i) = dsn_temp;
                self.dtau(:,i) = dtau_temp;
            end
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

        function GF = initialize_greens_functions(self, params, y, dx, varies_with_depth, variable_dip)
            [dx, y] = self.numerical_correction(dx, y, params, 1e-9);              
            xeval = y/(tan(params.dip*pi/180)) + dx;
            if and(~varies_with_depth, ~variable_dip)
                GF{1} = GreensFunctions(y);
                 if params.width_FW > 0  % and add a criterium to check if the FW dP and dT are not 0
                GF{1} = GF{1}.green_FW(xeval, y, params.dip, params.thick, params.throw, params.width_FW, 0, 0 );
            end
            if params.width_HW > 0 % and add a criterium to check if the FW dP and dT are not 0
                GF{1} = GF{1}.green_HW(xeval, y, params.dip, params.thick, params.throw, params.width_HW, 0, 0 );
            end
            elseif and(varies_with_depth, ~variable_dip)
                % calculate the same GF once, but shift it along depth axis
                slice_thick = (y(1) - y(2));                        % depth slice thickness
                y2 = [y; y(1:end-1)+(y(end)-y(1))-slice_thick];     % pad with zeros
                xeval2 = y2/(tan(params.dip*pi/180)) + dx;    
                i_mid = ceil(length(y2)/2);
                slice_y = y2(i_mid);                                % depth slice mid y on fault
                slice_x = slice_y/(tan(params.dip*pi/180));         % depth slice mid x on fault
                slice_throw = 0;                                    % depth slice throw, set to 0. 
                greens_f = GreensFunctions(y2);                % initialize Green's functions
                if params.width_FW > 0  % and add a criterium to check if the FW dP and dT are not 0
                    greens_f = greens_f.green_FW(xeval2, y2, params.dip, slice_thick, slice_throw, params.width_FW, slice_x, slice_y);
                end
                if params.width_HW > 0 % and add a criterium to check if the FW dP and dT are not 0
                    greens_f = greens_f.green_HW(xeval2, y2, params.dip, slice_thick, slice_throw, params.width_HW, slice_x, slice_y);
                end
                GF = cell(length(y), 1);
                for j = 1 : length(y)
                    GF{j} = greens_f;
                    GF{j}.Gnorm_FW = GF{j}.Gnorm_FW(i_mid-j+1:2*i_mid-j);
                    GF{j}.Gnorm_HW = GF{j}.Gnorm_HW(i_mid-j+1:2*i_mid-j);
                    GF{j}.Gshear_FW = GF{j}.Gshear_FW(i_mid-j+1:2*i_mid-j);
                    GF{j}.Gshear_HW = GF{j}.Gshear_HW(i_mid-j+1:2*i_mid-j);
                    % here, remove the xx components if needed
                end
            else
                % calculate separate GF at each depth interval
                slice_thick = y(1) - y(2);                      % depth slice thickness
                slice_throw = 0; 
                GF = cell(length(y), 1);
                % depth slice throw, set to 0. 
                for j = 1 : length(y)
                    slice_y = y(j); % - 0.5*slice_thick;        % depth slice mid y on fault
                    slice_x = slice_y/(tan(params.dip*pi/180)); % depth slice mid x on fault
                    GF{j} = GreensFunctions(y);              % initialize Green's functions
                    GF{j} = GF{j}.green_FW(xeval, y, params.dip, slice_thick, slice_throw, params.width_FW, slice_x, slice_y);
                    GF{j} = GF{j}.green_HW(xeval, y, params.dip,  slice_thick, slice_throw, params.width_HW, slice_x, slice_y );
                end
            end

        end

        function [dx, y] = numerical_correction(~, dx, y, params, correction_value)
            % numerical_correction Corrects x and y in case these coincide
            % with any of the boundaries of the reservoir shapes. this
            % results in singularities and division by 0
            if dx == 0 || dx == params.width_FW || dx == -params.width_HW
                dx = dx + correction_value;
            end
            reservoir_boundaries = [params.top_FW_y, params.top_HW_y, params.base_HW_y, params.base_FW_y];
            if any(ismember(y, reservoir_boundaries))
                i_boundary = find(ismember(y, reservoir_boundaries));
                y(i_boundary) = y(i_boundary) + correction_value;
            end
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
            vary_P = 0;
            vary_T = 0;
            if length(unique(pressure.dp_FW)) > size(pressure.dp_FW,2) + 1
                vary_P = 1;
            elseif length(unique(pressure.dp_HW)) > size(pressure.dp_FW,2) + 1
                vary_P = 1;
            end
        end
        % TODO: functionality to return stress in xx, yy, xy
        % TODO: temperature

        function [dscu] = get_scu_change(self, mu, cohesion)
            % return Shear Capacity Utilization 
            % TODO: check if this works for depth-dependent mu
            dscu = self.dtau ./ (self.dsn .* mu + cohesion); 
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
