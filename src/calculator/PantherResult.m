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
        load_table table
    end

    methods
        function self = PantherResult(y)
            self.y = y;
        end

        function self = make_result_summary(self, analysis)
            warning('off');
            for i = 1 : length(self.stress)
                self.summary.reactivation(i) = self.slip{i}.reactivation;
                self.summary.reactivation_load_step(i) = self.slip{i}.reactivation_load_step;
                self.summary.nucleation(i) = self.slip{i}.nucleation;
                self.summary.nucleation_load_step(i) = self.slip{i}.nucleation_load_step;
                n_steps = linspace(1,length(analysis.load_table.time_steps),length(analysis.load_table.time_steps));
                if ~isnan(self.slip{i}.reactivation_load_step)
                    if strcmp(analysis.load_case,'P')
                        self.summary.reactivation_dp(i) = interp1(n_steps, analysis.load_table.P_steps, self.slip{i}.reactivation_load_step);
                        self.summary.reactivation_dT(i) = nan;
                    elseif strcmp(analysis.load_case,'T') 
                        self.summary.reactivation_dT(i) = interp1(n_steps, analysis.load_table.T_steps, self.slip{i}.reactivation_load_step);
                        self.summary.reactivation_dp(i) = nan;
                    elseif strcmp(analysis.load_case,'PT')
                        self.summary.reactivation_dT(i) = interp1(n_steps, analysis.load_table.T_steps, self.slip{i}.reactivation_load_step);
                        self.summary.reactivation_dp(i) = nan;
                    end
                else
                        self.summary.reactivation_dp(i) = nan;
                        self.summary.reactivation_dT(i) = nan;
                end
                if ~isnan(self.slip{i}.nucleation_load_step)
                    if strcmp(analysis.load_case,'P') 
                        self.summary.nucleation_dp(i) = interp1(n_steps, analysis.load_table.P_steps, self.slip{i}.nucleation_load_step);
                        self.summary.nucleation_dT(i) = nan;
                    elseif strcmp(analysis.load_case,'T') 
                        self.summary.nucleation_dT(i) = interp1(n_steps, analysis.load_table.T_steps, self.slip{i}.nucleation_load_step);
                        self.summary.nucleation_dp(i) = nan;
                    elseif strcmp(analysis.load_case,'PT') 
                        self.summary.nucleation_dT(i) = interp1(n_steps, analysis.load_table.T_steps, self.slip{i}.nucleation_load_step);
                        self.summary.nucleation_dp(i) = nan;
                    end
                else
                    self.summary.nucleation_dp(i) = nan;
                    self.summary.nucleation_dT(i) = nan;
                end
%                 [self.summary.max_cff_rate(i), self.summary.mid_cff_rate(i)]  = self.get_cff_rates(analysis, i); 
            end
            warning('on'); 
        end       

    end

end
