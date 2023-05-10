classdef PantherPressure 
    % PantherPressure Computes pressure profiles along fault
    % Loes Buijze 13 - 04 - 2023
    
    properties
        p0
        dp_FW
        dp_HW
        dp_fault
    end

    properties (Dependent)
        p
    end

    methods
        function self = PantherPressure(member, y, loads, load_case, diffusion, p_fault_mode)
            ini = InitialPressure(member, y, p_fault_mode);                    % instance containing initial pressure
            self.p0 = ini.p0;                           % initial fault pressure
            p_steps = loads.P_steps;
            time_steps = loads.time_steps;
            [next_to_FW, next_to_HW, ~] = is_adjacent_to_reservoir(y, member.thick, member.throw);
            dp_unit_FW = double(next_to_FW);       % unit dp in FW compartment 
            dp_unit_HW = double(next_to_HW);       % unit dp in HW compartment
            self.dp_FW = dp_unit_FW * p_steps' .* loads.P_factor_FW';    % (dp, time) array of pressures in the footwall
            self.dp_HW = dp_unit_HW * p_steps' .* loads.P_factor_HW';    % (dp, time) array of pressures in the hanging wall
            
            % compute pressure diffusion to the seal and base
            % TODO enable different diffusivity in the fault zone
            if diffusion
                y_top = member.top_FW_y;
                y_base = member.base_FW_y;
                self.dp_FW = calc_dp_diffusion(y, y_top, y_base, time_steps, self.dp_FW, member.hyd_diffusivity);
                y_top = member.top_HW_y;
                y_base = member.base_HW_y;
                self.dp_HW = calc_dp_diffusion(y, y_top, y_base, time_steps, self.dp_HW, member.hyd_diffusivity);
            end
            
            % set the fault pressure w.r.t. HW and FW pressure
            if strcmp(p_fault_mode, 'max')
                self.dp_fault = max(self.dp_FW, self.dp_HW);
            elseif strcmp(p_fault_mode, 'min')
                self.dp_fault = min(self.dp_FW, self.dp_HW);
            elseif strcmp(p_fault_mode, 'mean')
                self.dp_fault = mean([self.dp_FW, self.dp_HW],2);
            elseif strcmp(p_fault_mode, 'FW')
                self.dp_fault = self.dp_FW;
            elseif strcmp(p_fault_mode, 'HW')
                self.dp_fault = self.dp_HW;
            end

            % plot pressures
%             plot(self.dp_HW(:, end), y, self.dp_FW(:, end), y);
%             hold on
%             plot(self.dp_fault(:, [1,5,10,25,50]), y, 'LineStyle','--','Color','k');
%             plot(self.dp_fault(:, end)/50, y, 'LineStyle','--','Color','b');
            
            self.dp_fault = member.p_factor_fault * self.dp_fault;
        end

       function p = get.p(self)
            p = self.p0 + self.dp_fault;
        end


    end

end