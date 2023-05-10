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
                max_cff = 0;        % maximum stress rate
                mid_cff = 0;        % stress rate at mid reservoir depth
                % TODO find better metric for mean stress
                cff = self.stress{i}.get_cff(self.ensemble.f_s(i), self.ensemble.cohesion(i));
                time = analysis.load_table.time_steps;      % time in yrs
                if self.summary.reactivation(i) && analysis.aseismic_slip
                    if self.summary.reactivation_index(i) > 1
                        cff = cff(:, 1:self.summary.reactivation_index(i));
                        time = time(1:self.summary.reactivation_index(i));
                        %                         tops = [analysis.ensemble{i}.top_FW_i(self.y), analysis.ensemble{i}.top_HW_i(self.y)];
%                         top = max(tops);
%                         bases = [analysis.ensemble{i}.base_FW_i(self.y), analysis.ensemble{i}.base_HW_i(self.y)];
%                         base = max(bases);
%                         mid_cff = mean(mean(cff_rate(top:base,:), 2));
                    end
                end
                cff_rate = diff(cff, [], 2) ./ diff(time)';
                max_cff = max(max(cff_rate));   % maximum Coulomb stress rate
                i_mid = ceil(length(self.y)/2);
                mid_cff = mean(cff_rate(i_mid,:));
           
        end

    end

end
