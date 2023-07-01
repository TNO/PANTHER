classdef Temperature
    % Sets the temperature changes in the hanging wall and footwall
    % compartments. 
    % T0:       initial temperature
    % dT_HW     temperature change hanging wall (len(y) x len(time))
    % dT_FW     temperature change footwall (len(y) x len(time))
    % dT_fault  temperature at fault (irrelevant for stress change)
    % Loes Buijze 13 - 04 - 2023
    
    properties
        T0
        dT_HW
        dT_FW
        dT_fault
    end

    properties (Dependent)
        T
    end

    methods
        function self = Temperature(member, y, loads, diffusion, T_fault_mode)
            % Initialize temperature loads
            self.T0 = -(y + member.depth_mid)*member.T_grad/1000 + member.T_offset;                           
            T_steps = loads.T_steps;               % degrees
            time_steps = loads.time_steps;         % years
            [next_to_FW, next_to_HW, ~] = is_adjacent_to_reservoir(y, member.thick, member.throw);
            dT_unit_FW = double(next_to_FW);       % unit dT in FW compartment 
            dT_unit_HW = double(next_to_HW);       % unit dT in HW compartment
            self.dT_FW = dT_unit_FW * T_steps' .* loads.T_factor_FW';    % (dp, time) array of pressures in the footwall
            self.dT_HW = dT_unit_HW * T_steps' .* loads.T_factor_HW';    % (dp, time) array of pressures in the hanging wall
            
            % compute temperature diffusion
            if diffusion
                y_top = member.top_FW_y;
                y_base = member.base_FW_y;
                self.dT_FW = calc_dT_diffusion(y, y_top, y_base, time_steps, self.dT_FW, member.therm_diffusivity);
                y_top = member.top_HW_y;
                y_base = member.base_HW_y;
                self.dT_HW = calc_dT_diffusion(y, y_top, y_base, time_steps, self.dT_HW, member.therm_diffusivity);
            end
            
            % set the fault temperature w.r.t. HW and FW temperature
            % irrelevant for stress changes. for plotting purposes. 
            if strcmp(T_fault_mode, 'max')
                self.dT_fault = max(self.dT_FW, self.dT_HW);
            elseif strcmp(T_fault_mode, 'min')
                self.dT_fault = min(self.dT_FW, self.dT_HW);
            elseif strcmp(T_fault_mode, 'mean')
                self.dT_fault = mean([self.dT_FW, self.dT_HW],2);
            elseif strcmp(T_fault_mode, 'FW')
                self.dT_fault = self.dT_FW;
            elseif strcmp(T_fault_mode, 'HW')
                self.dT_fault = self.dT_HW;
            end
            
        end

       function T = get.T(self)
            T = self.T0 + self.dT_fault;
        end


    end

end