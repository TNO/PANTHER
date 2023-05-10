classdef FaultStress
    % Fault stress object. Calculates and stores 
    % shear and normal stresses from initial and stress changes 

    properties
        sne double % effective normal stress (len(y), len(timesteps))
        tau double % shear stress (len(y), len(timesteps)). positive: normal faulting
    end

    properties (Dependent)
    end

    methods
        function self = FaultStress(size_y, size_t)
            % Initializes stress arrays
            self.sne = zeros(size_y, size_t);
            self.tau = zeros(size_y, size_t);
        end

        function self = compute_fault_stress(self, initial, change, pressure)
            self.sne = initial.sn0 + change.dsn - pressure;
            self.tau = initial.tau0 + change.dtau;
        end

        function scu = get_scu(self, mu, cohesion)
            % calculate the Shear Capacity Utilization
            scu = self.tau ./ (self.sne .* mu + cohesion); 
        end

        function cff = get_cff(self, mu, cohesion)
            % calculate the Coulomb Failure Stress
            cff = self.tau - (self.sne .* mu + cohesion); 
        end

        function [sne_f, tau_f] = get_failure_stress(self, i_fail)
            % obtain the fault stresses at arbitrary value between 1 and
            % length(timesteps)
            sne_f = zeros(size(self.sne, 1), 1);
            tau_f = zeros(size(self.sne, 1), 1);
            if i_fail < 1 || i_fail > size(self.sne, 2)
                disp('Selected time is outside calculation time range');
            else
                for i = 1 : size(self.sne, 1)
                    x_ind = linspace(1, size(self.sne, 2), size(self.sne, 2));% indices of time, P, T steps
                    sne_f(i) = interp1(x_ind, self.sne(i, :), i_fail);
                    tau_f(i) = interp1(x_ind, self.tau(i, :), i_fail);
                end     
            end
         end

    end
end