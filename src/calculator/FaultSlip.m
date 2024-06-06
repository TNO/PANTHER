classdef FaultSlip
    % Fault slip object. Contains shear slip, length of the fault that is
    % slipping, and identified when fault reactivation and nucleation occur

    properties
        slip double                 % shear slip on the fault
        slip_length double          % length of fault that is slipping (num rows is number of slip patches)
        reactivation logical        % indicator whether fault has been reactivated
        reactivation_load_step double   % load index at which fault has been reactivated
        nucleation logical          % indicator whether nucleation has occurred on the fault
        nucleation_load_step double     % load index at which nucleation occurs
    end

    properties (Dependent)
    end

    methods
        function self = FaultSlip(n_rows, n_cols)
            % initialize FaultSlip class.  
            self.slip = zeros(n_rows, n_cols);
            self.reactivation = 0;
            self.reactivation_load_step = nan(1,1);
            self.nucleation = 0;
            self.nucleation_load_step = nan(1,1);
        end

        function [self, tau_slip] = calculate_fault_slip(self, y, sne, tau, tau_f, mu_II)
            % calculate shear slip on fault
            tau_slip = zeros(size(sne));        % shear stress including stress redistribution due to slip
            K = self.stiffness_matrix(y, mu_II);
            for i = 1 : size(tau, 2)
                [tau_slip(:,i), islip{i}, slip_timestep] = calc_aseismic_slip(y, tau(:,i), tau_f(:,i), mu_II, K);
                self.slip(:,i) = slip_timestep;
            end
        end

        function self = reduce_steps(self, steps)
            props = properties(self);
            % reduce slip matrix and slip patch array
            for i = 1 : 2
                self.(props{i}) = self.(props{i})(:, steps);
            end             
        end

        function [self] = detect_nucleation(self, cell_length, sne, tau, f_s, f_d, d_c, cohesion, mu_II, nuc_crit, nuc_len_fixed)
            % identify whether nucleation occurs by comparing slip length
            % to theoretical nucleation length by Uenishi & Rice 2003
            % INPUT
            % sne   array size(y) x size(t). effective normal stress
            % tau   array size(y) x size(t). shear stress
            % f_s   static friction coefficient
            % f_d   dynamic friction coefficient
            % cohesion  cohesion
            % mu_II shear modulus mode II
            tau_f = sne.* f_s + cohesion;
            y2L = cell_length;                  % replace with dip
            slipping = (tau >= tau_f);
            % slip_zone_indices and ..length will have size(# slip zone x
            % length(time). 
            slip_zone_indices = self.get_slip_zone_indices(slipping);
            slip_zone_length = nan(size(slip_zone_indices));        % along-fault length of slip zones
            nucleation_length = nan(size(slip_zone_indices));       % nucleation length, per slip zone 
            for i = 1 : size(slip_zone_indices,2)
                for j = 1 : size(slip_zone_indices, 1)
                    if ~isempty(slip_zone_indices{j, i})
                        slip_zone_length(j,i) = (length(slip_zone_indices{j,i}) + 1) * y2L;
                        sne_slip = (sne(slip_zone_indices{j,i}, i));
                        tau_slip = (tau(slip_zone_indices{j,i}, i));
                        if length(f_s) == size(tau,1)
                            % add ,i here if friction changes per timestep
                            f_s_slip = f_s(slip_zone_indices{j,i}) ;
                        else
                            f_s_slip = f_s;
                        end
                        if length(f_d) == size(tau,1)
                            f_d_slip = (f_d(slip_zone_indices{j,i}));
                        else
                            f_d_slip = f_d;
                        end
                        if length(d_c) == size(tau,1)
                            d_c_slip = mean(d_c(slip_zone_indices{j,i}));
                        else
                            d_c_slip = d_c;
                        end
                        delta_tau = mean(sne_slip .* (f_s_slip - f_d_slip));
                        tau_0_d = mean(tau_slip - sne_slip .* f_d_slip);
                        nucleation_length(j,i) = self.calculate_nucleation_length( delta_tau, tau_0_d, d_c_slip, mu_II, nuc_crit, nuc_len_fixed);
                        %nucleation_length(j,i) = 1.158 * mu_II * d_c./((f_s - f_d) * average_sne_in_slip_zone);
                    end
                end
            end
            % detect reactivation and nucleation, per slip zone
            reactivation_per_slip_zone = repmat(size(slip_zone_indices,2), size(slip_zone_indices,1), 1);
            nucleation_per_slip_zone = nan(size(reactivation_per_slip_zone));  
            for j = 1 : size(slip_zone_indices, 1)
                if any(~isnan(nucleation_length(j,:))) && ~isempty(slip_zone_indices)
                    self.reactivation = 1;
                    reactivation_per_slip_zone(j) = find(~isnan(nucleation_length(j,:)),1,'first');
                    nucleation_length_j = nucleation_length(j, ~isnan(nucleation_length(j,:)));
                    slip_length_j = slip_zone_length(j, ~isnan(nucleation_length(j,:)));
                    l_difference = nucleation_length_j - slip_length_j;
                    x = linspace(1, length(l_difference), length(l_difference)) + reactivation_per_slip_zone(j) - 1;
                    if l_difference(1) < 0
                        % if the first timestep of reactivation already
                        % reaches nucleation, take that step
                        nucleation_per_slip_zone(j) = x(1);
                    else
                        if length(x) > 1 && ~strcmp(nuc_crit, 'fixed')
                            nucleation_per_slip_zone(j) = interp1(l_difference, x, 0);
                        elseif length(x) > 1 && strcmp(nuc_crit, 'fixed')
                            [~, unique_index] = unique(l_difference);
                            nucleation_per_slip_zone(j) = interp1(l_difference(unique_index), x(unique_index), 0);
                        end
                    end
                end
            end
            self.reactivation_load_step = min(reactivation_per_slip_zone);
            self.nucleation_load_step = min(nucleation_per_slip_zone);
            if ~isnan(self.nucleation_load_step)
                self.nucleation = 1;
            end
            self.slip_length = slip_zone_length;
        end

        function slip_zone_indices = get_slip_zone_indices(~, slipping)
            % identify different slip zones on the fault
            number_of_slip_zones = zeros(1,size(slipping,2));
            for i = 1 : size(slipping, 2)
                i_slip = find(slipping(:,i));
                if ~isempty(i_slip)
                    iz = find(diff(i_slip) > 1);            % gap in indices of slipping cells (multiple slip patches)
                    if size(iz,1) > 1
                        iz = iz';
                    end
                    if isempty(iz)
                        number_of_slip_zones(i) = 1;                               % number of slip zones
                        slip_zone_indices{1,i} = i_slip;
                    else
                        number_of_slip_zones(i) = length(iz)+1;                    % nuber of slip zones      
                        iz = [1, iz, length(i_slip)];                               % start and end indices of slip zones
                        for j = 1 : number_of_slip_zones(i)
                            if j == 1 
                                start_index = iz(1);
                            else
                                start_index = iz(j) + 1;
                            end
                            slip_zone_indices{j,i} = i_slip(start_index:iz(j+1)); 
                        end
                    end
                else
                     slip_zone_indices{1,i} = [];
                end
            end
        end

        function [K] = stiffness_matrix(~, y, mu_II)
            % compute stiffness matrix for aseismic slip (td_solve)
            nx = length(y);
            Kline = -nx/(2*pi*(y(1)-y(end))) ./ ( [0:nx-1]'.^2-0.25 );
            Ko = toeplitz(Kline);                   
            K = mu_II*Ko;
            K(and(K<0,K>(min(min(K)/10000)))) = 0;    % set very small changes to 0 to avoid continued interactions of the two peaks to 
        end


        function nuc_length = calculate_nucleation_length(~, delta_tau, tau_0_d, d_c, mu_II, nucleation_criterion, nucleation_length)
        % simplified comparison to theoretical nucleation lengths, based on
        % average values of sn, delta_tau within the slip zone
            if strcmp(nucleation_criterion, 'UR2D')
                nuc_length = 1.158 * mu_II * d_c./(delta_tau);  % delta_tau
            elseif strcmp(nucleation_criterion, 'Day3D')
                % Day 2005
                 nuc_length = 2* ( mu_II * d_c * ( delta_tau ) ) / ( pi * ( tau_0_d ) ^ 2 );
            elseif strcmp(nucleation_criterion, 'Ruan3D')
                % from script Vincent
                 nuc_length =  ((3.82 * pi)/4)^0.5 * (mu_II / (delta_tau))*d_c;
            elseif strcmp(nucleation_criterion, 'fixed')
                % fixed nucleation length
                 nuc_length =  nucleation_length;
            end
        end

%         function [d_max] = get_max_slip(self)
%         end
% 
%         function [d_nuc] = get_slip_at_nucleation(self)
%         end


    end
end
