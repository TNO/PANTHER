classdef FaultStress
    % Fault stress object. Calculates and stores 
    % shear and normal stresses from initial and stress changes 

    properties
        sne double      % effective normal stress (len(y), len(timesteps))
        tau double      % shear stress (len(y), len(timesteps)). positive: normal faulting
        sne_reac double % normal stress at reactivation (len(y), 1)
        tau_reac double % shear stress at reactivation (len(y), 1)
        sne_nuc double  % normal stress at nucleation (len(y), 1)
        tau_nuc double  % shear stress at nucleation (len(y), 1)
    end

    properties (Dependent)
    end

    methods
        function self = FaultStress(size_y, size_t)
            % Initializes stress arrays
            self.sne = zeros(size_y, size_t);
            self.tau = zeros(size_y, size_t);
            self.sne_reac = zeros(size_y, 1);
            self.tau_reac = zeros(size_y, 1);
            self.sne_nuc = zeros(size_y, 1);
            self.tau_nuc = zeros(size_y, 1);
        end

        function self = compute_fault_stress(self, initial, change, pressure)
            self.sne = initial.sn0 + change.dsn - pressure;
            self.tau = initial.tau0 + change.dtau;
        end


        function self = get_reactivation_nucleation_stress(self, reactivation_load_step, nucleation_load_step)
            [self.sne_reac, self.tau_reac] = self.get_stress_at_load_step(reactivation_load_step);
            [self.sne_nuc, self.tau_nuc] = self.get_stress_at_load_step(nucleation_load_step);
        end


        function scu = get_scu(self, f_s, cohesion)
            % calculate the Shear Capacity Utilization
            scu = self.tau ./ (self.sne .* f_s + cohesion); 
        end

        function cff = get_cff(self, mu, cohesion)
            % calculate the Coulomb Failure Stress
            cff = self.tau - (self.sne .* mu + cohesion); 
        end

        function [sne_f, tau_f] = get_stress_at_load_step(self, load_step)
            % obtain the fault stresses at arbitrary value between 1 and
            % length(timesteps)
            sne_f = zeros(size(self.sne, 1), 1);
            tau_f = zeros(size(self.sne, 1), 1);
            if load_step < 1 || load_step > size(self.sne, 2)
                disp('Selected time is outside calculation time range');
            else
                for i = 1 : size(self.sne, 1)
                    x_ind = linspace(1, size(self.sne, 2), size(self.sne, 2));% indices of time, P, T steps
                    sne_f(i) = interp1(x_ind, self.sne(i, :), load_step);
                    tau_f(i) = interp1(x_ind, self.tau(i, :), load_step);
                end     
            end
        end

        function [cff_max, cff_ymid] = get_cff_rates(self, f_s, cohesion, time_yrs, time_range, depth_range)
            % Get the maximum Coulomb Stress Change CFF rate along the
            % fault as well as the CFF
            % INPUT
            % time      % time in yrs. should be equal to size(self.sne, 2)
            if nargin < 4
                min_index  = time_range(1);
                max_index = time_range(2);    % must be indices within
            end
            cff_max = 0;        % maximum stress rate
            cff_ymid = 0;        % stress rate at mid reservoir depth
            if size(time_yrs, 2) ~= size(self.sne,2)
                error('Time steps do not match number of stress output times');
            end
            cff = self.get_cff(f_s, cohesion);
            cff = cff(:, min_index:max_index);
            time_yrs = time_yrs(min_index:max_index);
            cff_rate = diff(cff, [], 2) ./ diff(time_yrs)';
            cff_max = max(max(cff_rate));   % maximum Coulomb stress rate
            average_cff_over_depth_range = mean(cff_rate(depth_range,:));
            i_ymid = ceil(length(self.y)/2);
            cff_ymid = mean(cff_rate(i_ymid,:));
        end

        function self = reduce_steps(self, steps)
            props = properties(self);
            % iterate over non-dependent properties
            for i = 1 : 2
                self.(props{i}) = self.(props{i})(:, steps);
            end             
        end

    end
end