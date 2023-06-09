classdef FaultSlip
    % Fault slip object. Contains shear slip, length of the fault that is
    % slipping, and identified when fault reactivation and nucleation occur

    properties
        slip double                 % shear slip on the fault
        slip_length double          % length of fault that is slipping (num rows is number of slip patches)
        reactivation logical        % indicator whether fault has been reactivated
        reactivation_index double   % load index at which fault has been reactivated
        nucleation logical          % indicator whether nucleation has occurred on the fault
        nucleation_index double     % load index at which nucleation occurs
    end

    properties (Dependent)
    end

    methods
        function self = FaultSlip(n_rows, n_cols)
            % initialize FaultSlip class.  
            self.slip = zeros(n_rows, n_cols);
            self.reactivation = 0;
            self.reactivation_index = nan(1,1);
            self.nucleation = 0;
            self.nucleation_index = nan(1,1);
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

        function [self] = detect_nucleation(self, cell_length, sne, tau, f_s, f_d, d_c, cohesion, mu_II )
            % identify whether nucleation occurs by comparing slip length
            % to theoretical nucleation length by Uenishi & Rice 2003
            tau_f = sne.* f_s + cohesion;
            y2L = cell_length;% replace with dip
            slipping = (tau >= tau_f);
            slip_zone_indices = self.get_slip_zone_indices(slipping);
            slip_zone_length = nan(size(slip_zone_indices));        % along-fault length of slip zones
            nucleation_length = nan(size(slip_zone_indices));       % nucleation length, per slip zone 
            for i = 1 : size(slip_zone_indices,2)
                for j = 1 : size(slip_zone_indices, 1)
                    if ~isempty(slip_zone_indices{j, i})
                        slip_zone_length(j,i) = (length(slip_zone_indices{j,i}) + 1) * y2L;
                        nucleation_length(j,i) = 1.158 * mu_II * d_c./((f_s - f_d) * mean(sne(slip_zone_indices{j,i}, i)));
                    end
                end
            end
            % detect reactivation and nucleation, per slip zone
            reactivation_slip_zone = repmat(size(slip_zone_indices,2), size(slip_zone_indices,1), 1);
            nucleation_slip_zone = nan(size(reactivation_slip_zone));  
            for j = 1 : size(slip_zone_indices, 1)
                if any(~isnan(nucleation_length(j,:))) && ~isempty(slip_zone_indices)
                    self.reactivation = 1;
                    reactivation_slip_zone(j) = find(~isnan(nucleation_length(j,:)),1,'first');
                    nucleation_length_j = nucleation_length(j, ~isnan(nucleation_length(j,:)));
                    slip_length_j = slip_zone_length(j, ~isnan(nucleation_length(j,:)));
                    l_difference = nucleation_length_j - slip_length_j;
                    x = linspace(1, length(l_difference), length(l_difference)) + reactivation_slip_zone(j) -1;
                    if length(x) > 1
                        nucleation_slip_zone(j) = interp1(l_difference, x, 0);
                    end
                end
            end
            self.reactivation_index = min(reactivation_slip_zone);
            self.nucleation_index = min(nucleation_slip_zone);
            if ~isnan(self.nucleation_index)
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
                    iz = find(diff(i_slip) > 1);            % discontinuity in indices of slipping cells
                    if size(iz,1) > 1
                        iz = iz';
                    end
                    if isempty(iz)
                        number_of_slip_zones(i) = 1;                               % number of slip zones
                        slip_zone_indices{1,i} = i_slip;
                    else
                        number_of_slip_zones(i) = length(iz)+1;                    % nuber of slip zones      
                        iz = [1, iz, length(i_slip)];  % start and end indices of slipping zones
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

%         function [d_max] = get_max_slip(self)
%         end
% 
%         function [d_nuc] = get_slip_at_nucleation(self)
%         end


    end
end
