classdef PantherPressure
    
    % Loes Buijze 13 - 04 - 2023
    
    properties
        p0
        dp_fault
        thick
        throw
        diffusivity
        time_steps
        p_steps_reservoir
        p_fault_mode {mustBeMember(p_fault_mode,{'max','min','mean','FW','HW'})} = 'max';
        dp_fault_mode {mustBeMember(dp_fault_mode,{'max','max_abs','min', 'min_abs','mean','FW','HW'})} = 'min'; 
        p_res_mode {mustBeMember(p_res_mode, {'same','different'})} = 'same'
        diffusion_P logical = true 
    end

    properties (Dependent)
        p
    end

    methods
        % function self = PantherPressure(member, y, loads, load_case, diffusion, p_fault_mode, dp_fault_mode, p_res_mode)
        function self = PantherPressure(y, thick, throw, diffusivity, load_table, load_case, diffusion, p_fault_mode, dp_fault_mode, p_res_mode)
            member = struct();
            member.thick = thick; 
            member.throw = throw;
            self.thick = thick;
            self.throw = throw;
            self.diffusivity = diffusivity;
            self.time_steps = load_table.time_steps;
            self.p_steps_reservoir = load_table.P_steps;

            ini = InitialPressure(member, y, p_fault_mode, p_res_mode);                    % instance containing initial pressure
            self.p0 = ini.p0;                           % initial fault pressure
            self.p_steps_reservoir = loads.P_steps;
            time_steps = loads.time_steps;
            if contains(load_case, 'P')
                % if reservoir compartment on either side has 0 width, set
                % pressures in that compartment to 0. 
                if member.width_FW == 0
                    dp_FW = zeros(length(y), length(time_steps));       % set dP in FW compartment to 0
                else
                    p_side = 'FW';
                    dp_FW = self.get_pressure_change_on_side(y, member, ini.p0_FW, p_steps', loads.P_factor_FW', time_steps, diffusion, p_side);
                end
                if member.width_HW == 0
                    dp_HW = zeros(length(y), length(time_steps));       % set dP in HW compartment to 0
                else
                    p_side = 'HW';
                    dp_HW = self.get_pressure_change_on_side(y, member, ini.p0_HW, p_steps', loads.P_factor_HW', time_steps, diffusion, p_side);
                end
                self.dp_fault = zeros(size(dp_FW));
            else
                self.dp_fault = zeros(length(y), length(time_steps));
            end

            % set the fault depletion pressure w.r.t. HW and FW pressure
            if strcmp(dp_fault_mode, 'max')
                self.dp_fault = max(dp_FW, dp_HW);
            elseif strcmp(dp_fault_mode, 'max_abs')
                for i = 1 : size(dp_FW,2)
                    [~, ind] = max([abs(dp_FW(:,i)), abs(dp_HW(:,i))],[],2);
                    self.dp_fault(ind == 1, i) = dp_FW(ind == 1, i);
                    self.dp_fault(ind == 2, i) = dp_HW(ind == 2, i);
                end
            elseif strcmp(dp_fault_mode, 'min')
                self.dp_fault = min(dp_FW, dp_HW);
            elseif strcmp(dp_fault_mode, 'min_abs')
                for i = 1 : size(dp_FW,2)
                    [~, ind] = min([abs(dp_FW(:,i)), abs(dp_HW(:,i))],[],2);
                    self.dp_fault(ind == 1, i) = dp_FW(ind == 1, i);
                    self.dp_fault(ind == 2, i) = dp_HW(ind == 2, i);
                end
            elseif strcmp(dp_fault_mode, 'mean')
                self.dp_fault = mean([dp_FW, dp_HW],2);
            elseif strcmp(dp_fault_mode, 'FW')
                self.dp_fault = dp_FW;
            elseif strcmp(dp_fault_mode, 'HW')
                self.dp_fault = dp_HW;
            end
            
            self.dp_fault = member.p_factor_fault * self.dp_fault;
        end

        function [dp_on_side] = get_pressure_change_on_side(self, y, params, p0, p_steps, p_factor, time_steps, diffusion, p_side)
            p_on_side = self.get_pressure_on_side(y, params, p0, p_steps, p_factor, time_steps, diffusion, p_side);
            dp_on_side = p_on_side - p0;
        end

        function [p_on_side] = get_pressure_on_side(~, y, params, p0, p_steps, p_factor, time_steps, diffusion, p_side)
            thick = params.thick;
            throw = params.throw;
            if strcmp(p_side, 'FW')
                [next_to_p_side, ~, ~] = is_adjacent_to_reservoir(y, thick, throw);
                % y_base = -(thick - throw)/2;
                % y_top = (thick + throw)/2;
            elseif strcmp(p_side, 'HW')
                [~, next_to_p_side, ~] = is_adjacent_to_reservoir(y, thick, throw);
               %y_base = -(thick + throw)/2;
               % y_top = (thick - throw)/2; 
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end
            dp_unit = double(next_to_p_side);       % unit dp in FW compartment 
             if ~iscolumn(dp_unit)
                 dp_unit = dp_unit';
             end
             if ~isrow(p_steps)
                 p_steps = p_steps';
             end
             if ~isrow(p_factor)
                 p_factor = p_factor';
             end
             % dp_on_side = zeros(length(y), length(p_steps));
             dp_on_side = dp_unit * p_steps .* p_factor;    % (dp, time) array of pressures in the footwall
             p_on_side = dp_on_side + p0;
             if diffusion
                p_on_side = calc_dp_diffusion(y, y_top, y_base, time_steps, p_on_side, params.hyd_diffusivity);
             end
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

        function [p_at_load_step] = get_p_at_load_step(self, load_step, p_type)
            % obtain the pressure at certain loadstep, where load_step is
            % between 1 and height of loadtable
            self.check_if_property_exists(p_type);
            p_at_load_step = zeros(size(self.(p_type), 1), 1);
            if nargin < 3 
                p_type = 'p'; 
            end
            if load_step < 1 || load_step > size(self.(p_type), 2)
                disp('Selected load step is outside calculation time range');
            else
                for i = 1 : size(self.(p_type), 1)
                    x_ind = linspace(1, size(self.(p_type), 2), size(self.(p_type), 2));% indices of time, P, or T steps
                    p_at_load_step(i) = interp1(x_ind, self.(p_type)(i, :), load_step);
                end     
            end
        end

        function check_if_property_exists(self, queried_property_name)
            if ~ismember(queried_property_name, properties(self))
                error(['Queried property name is not a property of Class PantherPressure. ',newline,...
                    'Properties must be one of: ', strjoin(properties(self),', ')]);
            end

        end


        function plot(self)
            
            subplot(1,2,1)
            plot(self.p(:,1));
            hold on
            plot(self.p(:,end));
            subplot(1,2,2)
            plot(self.dp_fault(:,end));
     
        end


    end

end