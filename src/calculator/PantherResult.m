classdef PantherResult
    % PantherResult contains Panther output results

    properties
        y 
        pressure cell
        temperature cell
        stress cell
        slip cell
        summary table
        ensemble table
    end

    methods
        function self = PantherResult(y)
            self.y = y;
        end

        function self = make_result_summary(self, analysis)
            warning('off');
            for i = 1 : length(self.stress)
                self.summary.reactivation(i) = self.slip{i}.reactivation;
                self.summary.reactivation_index(i) = self.slip{i}.reactivation_index;
                self.summary.nucleation(i) = self.slip{i}.nucleation;
                self.summary.nucleation_index(i) = self.slip{i}.nucleation_index;
                if strcmp(analysis.load_case,'P') && ~isnan(self.slip{i}.nucleation_index)
                    n_steps = linspace(1,length(analysis.load_table.time_steps),length(analysis.load_table.time_steps));
                    self.summary.reactivation_dp(i) = interp1(n_steps, analysis.load_table.P_steps, self.slip{i}.reactivation_index);
                    self.summary.nucleation_dp(i) = interp1(n_steps, analysis.load_table.P_steps, self.slip{i}.nucleation_index);
                else
                    self.summary.reactivation_dp(i) = nan;
                    self.summary.nucleation_dp(i) = nan;
                end
                [self.summary.max_cff_rate(i), self.summary.mid_cff_rate(i)]  = self.get_cff_rates(analysis, i); 
            end
            warning('on'); 
        end

        function [max_cff, mid_cff] = get_cff_rates(self, analysis, i)
            % Get the maximum Coulomb Stress Change CFF rate along the
            % fault as well as the CFF at mid reservoir depth (y=0), 
            % averaged over the timesteps 
                max_cff = 0;        % maximum stress rate
                mid_cff = 0;        % stress rate at mid reservoir depth
                % TODO find better metric for mean stress
                cff = self.stress{i}.get_cff(self.ensemble.f_s(i), self.ensemble.cohesion(i));
                time = analysis.load_table.time_steps;      % time in yrs
                if self.summary.reactivation(i) && analysis.aseismic_slip
                    if self.summary.reactivation_index(i) > 1
                        cff = cff(:, 1:self.summary.reactivation_index(i));
                        time = time(1:self.summary.reactivation_index(i));
                    end
                end
                cff_rate = diff(cff, [], 2) ./ diff(time)';
                max_cff = max(max(cff_rate));   % maximum Coulomb stress rate
                i_mid = ceil(length(self.y)/2);
                mid_cff = mean(cff_rate(i_mid,:));
           
        end

        function [p_init] = get_initial_pressure(self)
            % Returns presssure at first timestep
            p_init = cell(length(self.pressure), 1);
            for i = 1 : length(self.pressure)
                p_init{i} = self.pressure{i}.p0(:,1);
            end
        end

        function [sne0, tau0] = get_initial_stress(self)
            % Returns shear and normal stress at first timestep
            sne0 = cell(length(self.stress), 1);
            tau0 = cell(length(self.stress), 1);
            for i = 1 : length(self.stress)
                sne0{i} = self.stress{i}.sne(:,1);
                tau0{i} = self.stress{i}.tau(:,1);
            end
        end

        function [dsne, dtau] = get_stress_changes(self)
            % Returns shear and normal stress change w.r.t. first time step
            dsne = cell(length(self.stress), 1);
            dtau = cell(length(self.stress), 1);
            for i = 1 : length(self.stress)
                dsne{i} = self.stress{i}.sne - self.stress{i}.sne(:,1);
                dtau{i} = self.stress{i}.tau - self.stress{i}.tau(:,1);
            end
        end


    end

end
