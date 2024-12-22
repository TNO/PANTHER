classdef (HandleCompatible) PantherPressure < ModelGeometry & FaultMesh
    
    
    properties
        hyd_diffusivity double = 2e-6
        time_steps (:,1) double = [0,1]     % [year] time steps 
        P_steps (:,1) = [0,-1]              % [MPa] pressure change
        P_factor_HW (:,1) = [1,1]
        P_factor_FW (:,1) = [1,1]
        P_factor_fault (:,1) = [1,1]
        p_fault_mode {mustBeMember(p_fault_mode,{'max','min','mean','FW','HW'})} = 'max';
        dp_fault_mode {mustBeMember(dp_fault_mode,{'max','max_abs','min', 'min_abs','mean','FW','HW'})} = 'min'; 
        p_res_mode {mustBeMember(p_res_mode, {'same','different'})} = 'same'
        diffusion_P logical = false
        p_grad (:,1) double {mustBePositive} = 10.2
        p_grad_res (:,1) double {mustBePositive} = 10.2
        p_offset (:,1) double  = 0
        p_over (:,1) double  = 0
    end

    properties (Dependent)
        p0
        p
        dp_fault
    end

    methods
        function self = PantherPressure(input_params, load_table, pressure_settings)
            if nargin >= 1
                self = self.update_properties(input_params);
            end
            if nargin >= 2
                self = self.update_properties(load_table);
            end
            if nargin >= 3
                self = self.update_properties(pressure_settings);
            end
        end  

       
        function [p0, p0_HW, p0_FW] = get_initial_pressure(self)
            [next_to_FW, next_to_HW, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw); 
            yy = self.y + self.depth_mid;    
            p0 = zeros(size(yy));                               % intialize fault pressure
            p0 = -(yy/1000).*self.p_grad + self.p_offset;       % [MPa] hydrostatic pressure
            p0_FW = p0; p0_HW = p0;                             % initialize FW and HW pressure with hydrostatic
            top_HW_i = self.top_HW_i(self.y);                      % index where HW compartment starts (top)
            top_FW_i = self.top_FW_i(self.y);                      % index where FW compartment starts (top)
            top_res_i = min(top_HW_i, top_FW_i);                % top most depth of reservoir interval
            base_FW_i = self.base_FW_i(self.y); 
            base_HW_i = self.base_HW_i(self.y);
            % set overpressure w.r.t. hydrostatic gradient, in reservoir
            % and base
            p0_HW(top_HW_i:end) = p0_HW(top_HW_i:end) + self.p_over;  
            p0_FW(top_FW_i:end) = p0_FW(top_FW_i:end) + self.p_over;  
            % set reservoir pressure gradient within the reservoir compartments
            if strcmp(self.p_res_mode, 'same')
                base_res_i = max(base_FW_i, base_HW_i );
                p0_FW(top_FW_i:base_res_i) = p0_FW(base_res_i) - (yy(top_FW_i:base_res_i) - yy(base_res_i))*self.p_grad_res/1000;
                p0_HW(top_HW_i:base_res_i) = p0_HW(base_res_i) - (yy(top_HW_i:base_res_i) - yy(base_res_i))*self.p_grad_res/1000;
            else
                p0_FW(next_to_FW) = p0_FW(base_FW_i) - (yy(next_to_FW) - yy(base_FW_i))*self.p_grad_res/1000;
                p0_HW(next_to_HW) = p0_HW(base_HW_i) - (yy(next_to_HW) - yy(base_HW_i))*self.p_grad_res/1000; 
            end
            p0 = self.set_fault_pressure(p0_HW, p0_FW, self.p_fault_mode);
        end

        function dp_fault = get_dp_fault(self) 

            dp_fault = zeros(length(self.y), length(self.time_steps));
            % if reservoir compartment on either side has 0 width, set
            % pressures in that compartment to 0. 
            if self.width_FW == 0
                dp_FW = zeros(length(self.y), length(self.time_steps));       % set dP in FW compartment to 0
            else
                % dp_FW = self.get_pressure_change_on_side(self.y, input_params, ini.p0_FW, p_steps', loads.P_factor_FW', time_steps, diffusion, p_side);
                dp_FW = self.get_pressure_change_on_side(self.P_factor_FW,'FW');
            end
            if self.width_HW == 0
                dp_HW = zeros(length(self.y), length(self.time_steps));       % set dP in HW compartment to 0
            else
                dp_HW = self.get_pressure_change_on_side(self.P_factor_HW,'HW');
            end
            dp_fault = self.set_fault_pressure(dp_HW, dp_FW, self.dp_fault_mode);
            dp_fault = self.P_factor_fault' .* dp_fault;
        end

        function [fault_pressure] = set_fault_pressure(self, p_HW, p_FW, side_mode)
        %set the fault depletion pressure w.r.t. HW and FW pressure
            fault_pressure = zeros(size(p_HW));
            if strcmp(side_mode, 'max')
                fault_pressure = max(p_FW, p_HW);
            elseif strcmp(side_mode, 'max_abs')
                for i = 1 : size(p_FW,2)
                    [~, ind] = max([abs(p_FW(:,i)), abs(p_HW(:,i))],[],2);
                    fault_pressure(ind == 1, i) = p_FW(ind == 1, i);
                    fault_pressure(ind == 2, i) = p_HW(ind == 2, i);
                end
            elseif strcmp(side_mode, 'min')
                fault_pressure = min(p_FW, p_HW);
            elseif strcmp(self.side_mode, 'min_abs')
                for i = 1 : size(p_FW,2)
                    [~, ind] = min([abs(p_FW(:,i)), abs(p_HW(:,i))],[],2);
                    fault_pressure(ind == 1, i) = p_FW(ind == 1, i);
                    fault_pressure(ind == 2, i) = p_HW(ind == 2, i);
                end
            elseif strcmp(side_mode, 'mean')
                for i = 1 : size(p_FW,2)
                    fault_pressure(:,i) = mean([p_FW(:,i), p_HW(:,i)],2);
                end
            elseif strcmp(side_mode, 'FW')
                fault_pressure = p_FW;
            elseif strcmp(side_mode, 'HW')
                fault_pressure = p_HW;
            end
        end

       
        function  [dp_on_side] = get_pressure_change_on_side(self, p_factor, p_side)
            % p_on_side = self.get_pressure_on_side(self.y, self, self.p0, self.P_steps, p_factor, self.time_steps, self.diffusion_P, p_side);
            p_on_side = self.get_pressure_on_side(p_factor, p_side);
            dp_on_side = p_on_side - self.p0;
        end

        function  [dp_HW] = get_HW_pressure_change(self)
            % p_on_side = self.get_pressure_on_side(self.y, self, self.p0, self.P_steps, p_factor, self.time_steps, self.diffusion_P, p_side);
            p_HW = self.get_HW_pressure();
            [~, p0_HW, ~] = self.get_initial_pressure();
            dp_HW = p_HW - p0_HW;    
        end


        function  [dp_FW] = get_FW_pressure_change(self)
            % p_on_side = self.get_pressure_on_side(self.y, self, self.p0, self.P_steps, p_factor, self.time_steps, self.diffusion_P, p_side);
            p_FW = self.get_FW_pressure();
            [~, ~, p0_FW] = self.get_initial_pressure();
            dp_FW = p_FW - p0_FW;    
        end

        function [p_HW] = get_HW_pressure(self)
            p_HW = self.get_pressure_on_side(self.P_factor_HW, 'HW');
        end

        function [p_FW] = get_FW_pressure(self)
            p_FW = self.get_pressure_on_side(self.P_factor_FW, 'FW');
        end


        function [p_on_side] = get_pressure_on_side(self, p_factor, p_side)
            % get_pressure_on_side Method to calculate pressure on a given side (FW or HW).
            %
            % This method calculates the pressure in a given reservoir compartment (FW or HW) based on
            % the input pressure factor for that compartment and side indicator. It also accounts for diffusion
            % if specified.
            %
            % Parameters:
            %   p_factor - Pressure factor for the side.
            %   p_side - Side indicator ('FW' or 'HW').
            %
            % Returns:
            %   p_on_side - Calculated pressure on the specified side.

            % Validate side indicator
            if strcmp(p_side, 'FW')
                [next_to_p_side, ~, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw);
            elseif strcmp(p_side, 'HW')
                [~, next_to_p_side, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw);
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end

            % Get top and base y-coordinates for the reservoir compartment
            y_top = self.get_top_y(p_side);
            y_base = self.get_base_y(p_side);

            % Ensure dp_unit is a column vector
            dp_unit = double(next_to_p_side);       
            if ~iscolumn(dp_unit)
                 dp_unit = dp_unit';
            end

            % Ensure P_steps and p_factor are row vectors
            p_factor = p_factor(:)';
             
             dp_on_side = dp_unit * self.P_steps' .* p_factor;
             p_on_side = dp_on_side + self.p0;
             if self.diffusion_P
                p_on_side = calc_dp_diffusion(self.y, y_top, y_base, self.time_steps, p_on_side, self.hyd_diffusivity);
             end
        end

       function a = get.p0(self)
            a =  self.get_initial_pressure();
       end
       
       function a = get.dp_fault(self)
            a = self.get_dp_fault();
       end

       function p = get.p(self)
            p = self.p0 + self.dp_fault;
       end

        function self = reduce_steps(self, steps)
            props = properties(self);
            % iterate over non-dependent properties
            if isnan(steps)     
                for i = 1 : length(props)
                    if ismember(props{i}, {'time_steps','P_steps','P_factor_HW', 'P_factor_FW','P_factor_fault'})
                        self.(props{i}) = [];
                    end             
                end
            else
                for i = 1 : length(props)
                    if ismember(props{i},{'time_steps','P_steps','P_factor_HW', 'P_factor_FW','P_factor_fault'})
                        self.(props{i}) = self.(props{i})(steps, :);
                    end             
                end
            end
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
                error(['Queried property name is not a property of Class Pressure. ',newline,...
                    'Properties must be one of: ', strjoin(properties(self),', ')]);
            end

        end

        function [self] = update_properties(self, input)
            if istable(input)
                self = self.updatePropertiesFromTable(input);
            elseif isobject(input)
                self = self.updatePropertiesFromClass(input);
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

        function self = updatePropertiesFromTable(self, inputTable)
            % Get the list of properties of the class
            props = properties(self);
            
            % Get the list of column names from the input table
            tableColumns = inputTable.Properties.VariableNames;
            
            % Loop through each property and update if there's a matching column in the table
            for i = 1:length(props)
                propName = props{i};
                if ismember(propName, tableColumns)
                    self.(propName) = inputTable.(propName);
                end
            end
        end

        function self = updatePropertiesFromClass(self, inputClass)
            % Get the list of properties of the class
            props = properties(self);
            
            % Get the list of column names from the input table
            property_names_of_input_class = properties(inputClass);
            
            % Loop through each property and update if there's a matching column in the table
            for i = 1:length(props)
                property_name = props{i};
                if ismember(property_name, property_names_of_input_class)
                    self.(property_name) = inputClass.(property_name);
                end
            end
        end

    end

end