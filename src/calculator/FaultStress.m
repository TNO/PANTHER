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
        end

        function self = get_nucleation_stress(self, nucleation_load_step)
            % INPUT
            % nucleation_load_step double 
            [self.sne_nuc, self.tau_nuc] = self.get_stress_at_load_step(nucleation_load_step);
        end

        function self = get_reactivation_stress(self, reactivation_load_step)
            [self.sne_reac, self.tau_reac] = self.get_stress_at_load_step(reactivation_load_step);
        end


        function [sne_f, tau_f] = get_stress_at_load_step(self, load_step)
            % obtain the fault stresses at arbitrary value between 1 and
            % length(timesteps). 
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

        function [dsne, dtau] = get_stress_changes(self)
            % Returns shear and normal stress change w.r.t. first time step
            dsne = self.sne - self.sne(:,1);
            dtau = self.tau - self.tau(:,1);
        end

        function self = reduce_steps(self, steps)
            props = properties(self);
            % iterate over non-dependent properties
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

    end
end