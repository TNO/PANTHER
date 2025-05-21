classdef FaultSlip
    % Fault slip object. Contains shear slip, length of the fault that is
    % slipping, and identified when fault reactivation and nucleation occur

    properties
        slip double                 % shear slip on the fault
        slip_length double          % length of fault that is slipping (num rows is number of slip patches)
        slip_zone_ymid double       % center depth of the slip zones
        reactivation logical        % indicator whether fault has been reactivated
        reactivation_load_step double   % load index at which fault has been reactivated
        nucleation logical          % indicator whether nucleation has occurred on the fault
        nucleation_load_step double % load index at which nucleation occurs
        nucleation_length double    % nucleation length at the nucleation point
        nucleation_zone_ymid double      % center y of the nucleation zone
        max_slip_length double
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
            self.nucleation_length = nan(1,1);
            self.nucleation_zone_ymid = nan(1,1);
            self.max_slip_length = nan(1,1);
        end

        function [self, tau_slip] = calculate_fault_slip(self, L, sne, tau, tau_f, mu_II)
            % calculate shear slip on fault
            % INPUT
            % L         [m] along-fault length
            % sne       [MPa] effective normal stress
            % tau       [MPa] shear stress
            % tauf_f    [MPa] fault shear strength (sne*f_s + c)    
            % mu_II     [MPa] mode II shear modulus
            % OUTPUT
            % tau_slip  
            tau_slip = zeros(size(sne));        % shear stress including stress redistribution due to slip
            % detect non-uniformly spaced along-fault length (which happens with varying
            % dip)
            non_uniform_L = 0;
            dL = diff(L);
            dL = round(dL*1000)/1000;      % remove tiny differences in spacing
            % for non-uniformly spaced L, set up new regularly spaced L
            if length(unique(dL)) > 1
                L_new = ones(size(L));
                L_new = L_new.*linspace(L(1), L(end), length(L))';
                non_uniform_L = 1;
            else
                L_new = L;
            end
            %L = self.
            K = self.stiffness_matrix(L_new, mu_II);
            islip = cell(size(tau,2),1);  % indices of slipping cells
            for i = 1 : size(tau, 2)
                % interpolate shear stress and strength for new L
                if non_uniform_L
                   tau_new = interp1(L,  tau(:,i), L_new);
                   tau_f_new = interp1(L, tau_f(:,i), L_new);
                else
                    tau_new = tau(:,i);
                    tau_f_new = tau_f(:,i);
                end
                [tau_slip_new, islip{i}, slip_timestep_new] = calc_aseismic_slip(L_new, tau_new, tau_f_new, mu_II, K);
                if non_uniform_L
                   tau_slip(:,i) = interp1(L_new, tau_slip_new, L);
                   slip_timestep = interp1(L_new, slip_timestep_new, L);
                else
                    slip_timestep = slip_timestep_new;
                    tau_slip(:,i) = tau_slip_new;
                end                
                self.slip(:,i) = slip_timestep;
            end
        end


        function self = reduce_steps(self, steps)
            props = properties(self);
            % reduce slip matrix and slip patch array
            if isnan(steps)
                for i = 1 : 2
                    self.(props{i}) = [];
                end 
            else
                for i = 1 : 2
                    self.(props{i}) = self.(props{i})(:, steps);
                end 
            end
        end

        function [self] = detect_nucleation(self, y, L, sne, tau, f_s, f_d, d_c, cohesion, mu_II, nuc_crit, nuc_len_fixed)
            % identify whether nucleation occurs by comparing slip length
            % to theoretical nucleation length by Uenishi & Rice 2003
            % INPUT
            % y     depth
            % L     along fault length L
            % sne   array size(y) x size(t). effective normal stress
            % tau   array size(y) x size(t). shear stress
            % f_s   static friction coefficient
            % f_d   dynamic friction coefficient
            % cohesion  cohesion
            % mu_II shear modulus mode II
            tau_f = sne.* f_s + cohesion;
            %y2L = cell_length;                  % replace with dip
            y2L = L;
            slipping = (tau >= tau_f);
            % slip_zone_indices and ..length will have size(number of slip zone x
            % length(time). 
            slip_zone_indices = self.get_slip_zone_indices(slipping);
            slip_zone_length = nan(size(slip_zone_indices));        % along-fault length of slip zones
            nucleation_length_per_slip_zone = nan(size(slip_zone_indices));       % nucleation length, per slip zone
            % iterate over time steps
            for i = 1 : size(slip_zone_indices,2)
                % iterate over number of slip zones
                for j = 1 : size(slip_zone_indices, 1)
                    if ~isempty(slip_zone_indices{j, i})
                        L_start = L(slip_zone_indices{j,i}(1));
                        L_end = L(slip_zone_indices{j,i}(end));
                        if L_start ~= L(1)
                            L_start = L_start - 0.5 * (L_start - L(slip_zone_indices{j,i}(1) - 1)); % add half element 
                        end
                        if L_end ~= L(end)
                            L_end = L_end + 0.5 * (L_start - L(slip_zone_indices{j,i}(1) - 1)); % add half element 
                        end
                        slip_zone_length(j,i) = L_start - L_end;
                        sne_slip = (sne(slip_zone_indices{j,i}, i));    % normal stresses at all slipping indices within a slip zone j
                        tau_slip = (tau(slip_zone_indices{j,i}, i));    % shear stresses at all within a slip zone j
                        
                        if length(f_s) == size(tau,1)
                            % add ,i here if friction changes per timestep
                            f_s_slip = f_s(slip_zone_indices{j,i}) ;    % friction at all slipping indices (heterogeneous f_s)
                        else
                            f_s_slip = f_s;                             % friction coefficient at all slipping indices (uniform f_s)
                        end
                        if length(f_d) == size(tau,1)
                            f_d_slip = (f_d(slip_zone_indices{j,i}));   % dynamic friction coefficient at all slipping indiced (heterogeneous f_d)
                        else
                            f_d_slip = f_d;                             % dynamic friction coefficient at all slipping indices (uniform f_d)
                        end
                        if length(d_c) == size(tau,1)
                            d_c_slip = mean(d_c(slip_zone_indices{j,i}));
                        else
                            d_c_slip = d_c;
                        end
                        if length(mu_II) == size(tau,1)
                            mu_II_slip = mean(mu_II(slip_zone_indices{j,i}));
                        else
                            mu_II_slip = mu_II;
                        end
                        delta_tau = mean(sne_slip .* (f_s_slip - f_d_slip));   % strength drop
                        tau_0_d = mean(tau_slip - sne_slip .* f_d_slip);        % mean potential stress drop in slip zone

                        nucleation_length_per_slip_zone(j,i) = self.calculate_nucleation_length(delta_tau, tau_0_d, d_c_slip, mu_II_slip, nuc_crit, nuc_len_fixed);
                        %nucleation_length(j,i) = 1.158 * mu_II * d_c./((f_s - f_d) * average_sne_in_slip_zone);
                    end
                end
            end
            % detect reactivation and nucleation, per slip zone
            % initialize reactivation and nucleation load step arrays 
            % (size n_slip zones x time_steps)
            reactivation_per_slip_zone = nan(size(slip_zone_indices,1), 1);
            nucleation_index_per_slip_zone = nan(size(reactivation_per_slip_zone));  
            y_mid = nan(size(slip_zone_indices));
            for j = 1 : size(slip_zone_indices, 1)
                if any(~isnan(nucleation_length_per_slip_zone(j,:))) && ~isempty(slip_zone_indices)
                    self.reactivation = 1;
                    reactivation_per_slip_zone(j) = find(~isnan(nucleation_length_per_slip_zone(j,:)),1,'first');
                    nucleation_length_j = nucleation_length_per_slip_zone(j, ~isnan(nucleation_length_per_slip_zone(j,:)));
                    slip_length_j = slip_zone_length(j, ~isnan(nucleation_length_per_slip_zone(j,:)));
                    l_difference = nucleation_length_j - slip_length_j;
                    x = linspace(1, length(l_difference), length(l_difference)) + reactivation_per_slip_zone(j) - 1;
                    if l_difference(1) < 0
                        % if the first timestep of reactivation already
                        % reaches nucleation, take that step
                        nucleation_index_per_slip_zone(j) = x(1);
                    else
                        if any(l_difference < 0)
                            [~, unique_index] = unique(l_difference);
                            nucleation_index_per_slip_zone(j) = interp1(l_difference(unique_index), x(unique_index), 0);
                        end
                    end
                end
                for i = 1 : size(slip_zone_indices,2) 
                    if ~isempty(slip_zone_indices{j, i})
                        y_tops(j, i) = y(slip_zone_indices{j, i}(1));
                        y_base(j, i) = y(slip_zone_indices{j, i}(end));
                        y_mid(j, i) = (y_tops(i) + y_base(i))/2;
                    end
                end
                self.slip_zone_ymid = y_mid;
                
            end
            self.slip_length = slip_zone_length;
            % find the earliest reactivation and nucleation for the
            % different slip zones
            self.reactivation_load_step = min(reactivation_per_slip_zone);
            [self.nucleation_load_step, nucleation_zone_number] = min(nucleation_index_per_slip_zone);
            if ~isnan(self.nucleation_load_step)
                self.nucleation = 1;
                indices = linspace(1, length(nucleation_length_per_slip_zone),length(nucleation_length_per_slip_zone));
                % get the nucleation length (alternatively derive from
                % slip_zone_length and nucleation_load_step)
                self.nucleation_length = interp1(indices, nucleation_length_per_slip_zone(nucleation_zone_number,:), self.nucleation_load_step);
                % nucleation_zone_indices = slip_zone_indices(nucleation_zone_number,:);
                self.nucleation_zone_ymid =  interp1(indices, self.slip_zone_ymid(nucleation_zone_number,:), self.nucleation_load_step);
                self.max_slip_length = self.nucleation_length;
            else
                self.nucleation_zone_ymid = nan; 
                if self.reactivation
                    self.max_slip_length = max(max(self.slip_length));
                else
                    self.max_slip_length = nan;
                end
            end
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

        function [K] = stiffness_matrix(~, L, mu_II)
            % set stiffness matrix for aseismic slip (td_solve)
            % INPUT
            % L         along fault length
            % nx        number of fault elements
            % mu_II     mode II shear stiffness
            nx = length(L);
            Kline = -nx/(2*pi*(L(1)-L(end))) ./ ( [0:nx-1]'.^2-0.25 );
            Ko = toeplitz(Kline);
            % K = Ko*mu_II;     % if stiffness not depth dependent
            K = Ko .* mu_II';   % take depth dependent stiffness into account
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

    end
end
