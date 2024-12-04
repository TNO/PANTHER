classdef PantherPressure
    
    % Loes Buijze 13 - 04 - 2023
    
    properties
        p0
        dp_fault
        dp_HW
        dp_FW
    end

    properties (Dependent)
        p
    end

    methods
        function self = PantherPressure(member, y, loads, load_case, diffusion, p_fault_mode, p_res_mode)
            ini = InitialPressure(member, y, p_fault_mode, p_res_mode);                    % instance containing initial pressure
            self.p0 = ini.p0;                           % initial fault pressure
            p_steps = loads.P_steps;
            time_steps = loads.time_steps;
            if contains(load_case, 'P')
                [next_to_FW, next_to_HW, ~] = is_adjacent_to_reservoir(y, member.thick, member.throw);
                dp_unit_FW = double(next_to_FW);       % unit dp in FW compartment 
                dp_unit_HW = double(next_to_HW);       % unit dp in HW compartment
                % if reservoir compartment on either side has 0 width, set
                % pressures in that compartment to 0. 
                if member.width_HW == 0
                    dp_unit_HW = zeros(size(dp_unit_HW));       % set dP in HW compartment to 0
                end
                if member.width_FW == 0
                    dp_unit_FW = zeros(size(dp_unit_FW));       % set dP in FW compartment to 0
                end
                self.dp_FW = dp_unit_FW * p_steps' .* loads.P_factor_FW';    % (dp, time) array of pressures in the footwall
                self.dp_HW = dp_unit_HW * p_steps' .* loads.P_factor_HW';    % (dp, time) array of pressures in the hanging wall
                
                p_FW_nodiffusion = self.dp_FW + ini.p0_FW;
                p_HW_nodiffusion = self.dp_HW + ini.p0_HW;
                % compute pressure diffusion to the seal and base
                if diffusion
                    % compute diffusion. employs difference in pressure
                    % between top and base FW and seal and base
                    y_top = member.top_FW_y;
                    y_base = member.base_FW_y;
                    p_FW_diffusion = calc_dp_diffusion(y, y_top, y_base, time_steps, p_FW_nodiffusion, member.hyd_diffusivity);
                    self.dp_FW =  p_FW_diffusion - ini.p0;
                    
                    % compute diffusion. employs difference in pressure
                    % between top and base HW and seal and base
                    y_top = member.top_HW_y;
                    y_base = member.base_HW_y;
                    p_HW_diffusion = calc_dp_diffusion(y, y_top, y_base, time_steps, p_HW_nodiffusion, member.hyd_diffusivity);
                    self.dp_HW =  p_HW_diffusion - ini.p0_HW;
                end

            else
                self.dp_HW = zeros(length(y), length(time_steps));
                self.dp_fault = zeros(length(y), length(time_steps));
                self.dp_FW = zeros(length(y), length(time_steps));
            end

            % set the fault depletion pressure w.r.t. HW and FW pressure
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
            
            self.dp_fault = member.p_factor_fault * self.dp_fault;
            
            h1 = figure(1); clf(h1);
            plot(self.dp_fault(:,:), y, self.p(:,:), y);

        end

        function self = reduce_steps(self, steps)
            props = properties(self);
            % iterate over non-dependent properties
            if isnan(steps)
                for i = 1 : length(props) - 1
                    if ~strcmp(props{i},'p0')
                        self.(props{i}) = [];
                    end             
                end
            else
                for i = 1 : length(props) - 1
                    if ~strcmp(props{i},'p0')
                        self.(props{i}) = self.(props{i})(:, steps);
                    end             
                end
            end
        end

       function p = get.p(self)
            p = self.p0 + self.dp_fault;
       end


    end

end