classdef (HandleCompatible) Pressure < ModelGeometry & FaultMesh
    % PantherPressure Class to represent the pressure calculations in the geological model.
    % A simplified pressure profile is prescribed, with uniform pressure
    % change in the reservoir and difusion to the seal and base (optional)
    % This class extends ModelGeometry and FaultMesh to include properties and methods
    % for calculating pressure changes in the model, including initial pressure, pressure
    % changes across faults, and pressure at specific load steps.
    %
    % Properties:
    %   hyd_diffusivity - Hydraulic diffusivity
    %   time_steps - Time steps for the simulation
    %   P_steps - Pressure change steps
    %   P_factor_HW - Pressure factor for the hanging wall
    %   P_factor_FW - Pressure factor for the footwall
    %   P_factor_fault - Pressure factor for the fault
    %   p_fault_mode - Mode for fault pressure calculation
    %   dp_fault_mode - Mode for fault pressure difference calculation
    %   p_res_mode - Mode for reservoir pressure calculation
    %   diffusion_P - Logical flag for diffusion
    %   p_grad - Pressure gradient
    %   p_grad_res - Reservoir pressure gradient
    %   p_offset - Pressure offset
    %   p_over - Overpressure
    %   load_case - Load case identifier
    %
    % Dependent Properties:
    %   p0 - Initial pressure
    %   p - Current pressure
    %   dp_fault - Pressure difference across the fault
    %
    % Methods:
    %   PantherPressure - Constructor to initialize the pressure model
    %   get_initial_pressure - Returns the initial pressure in the model
    %   get_dp_fault - Returns the pressure difference across the fault
    %   set_fault_pressure - Sets the fault pressure based on HW and FW pressures
    %   get_initial_pressure_on_side - Returns the HW or FW initial pressure 
    %   get_HW_p0 - Returns initial pressure in the hanging wall
    %   get_FW_p0 - Returns initial pressure in the footwall
    %   get_pressure_change_on_side - Calculates pressure change on a given side
    %   get_HW_pressure_change - Returns pressure change in the hanging wall
    %   get_FW_pressure_change - Returns pressure change in the footwall
    %   get_HW_pressure - Returns pressure in the hanging wall
    %   get_FW_pressure - Returns pressure in the footwall
    %   get_pressure_on_side - Calculates pressure on a given side
    %   reduce_steps - Reduces the number of steps in the simulation
    %   get_p_at_load_step - Returns pressure at a specific load step
    %   get_dp_at_load_step - Returns pressure difference at a specific load step
    %   check_if_property_exists - Checks if a property exists in the class
    %   update_properties - Updates properties from input parameters
    %   plot - Plots the pressure and pressure difference
    %   updatePropertiesFromTable - Updates properties from a table
    %   updatePropertiesFromClass - Updates properties from another class

    properties
        hyd_diffusivity double = 2e-6
        time_steps (:,1) double = [0,1]     % [year] time steps 
        P_steps (:,1) = [0,-1]              % [MPa] pressure change
        P_factor_HW (:,1) = [1,1]
        P_factor_FW (:,1) = [1,1]
        P_factor_fault (:,1) = [1,1]
        p0_fault_mode {mustBeMember(p0_fault_mode,{'max','min','mean','FW','HW'})} = 'max';
        p_fault_mode {mustBeMember(p_fault_mode,{'max','max_abs','min', 'min_abs','mean','FW','HW'})} = 'min'; 
        p_res_mode {mustBeMember(p_res_mode, {'same','different'})} = 'same'
        diffusion_P logical = false
        p_grad (:,1) double {mustBePositive} = 10.2
        p_grad_res (:,1) double {mustBePositive} = 10.2
        p_offset (:,1) double  = 0
        p_over (:,1) double  = 0
        load_case = 'P'
    end

    properties (Access=private)
    end

    properties (Dependent)
        p0
        p
        dp_fault
    end

    methods
        function self = Pressure(input_params, load_table, pressure_settings)
            % PantherPressure Constructor to initialize the pressure model.
            % Input:
            %   input_params - Input parameters for the model
            %   load_table - Load table for the model
            %   pressure_settings - Pressure settings for the model
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

       
        function [p0] = get_initial_pressure(self)
            % get_initial_pressure Returns the initial fault pressure
            % Output:
            %   p0 - Initial pressure in the fault
            p0_FW = self.get_FW_p0();
            p0_HW = self.get_HW_p0();
            p0 = self.set_fault_pressure(p0_HW, p0_FW, self.p0_fault_mode);
        end

       

        % function dp_fault = get_dp_fault(self) 
        %     % get_dp_fault Returns the pressure difference within the fault.
        %     % Output:
        %     %   dp_fault - Pressure difference within the fault
        %     % if reservoir compartment on either side has 0 width, set
        %     % pressures in that compartment to 0. 
        %     if self.width_FW == 0
        %         dp_FW = zeros(length(self.y), length(self.time_steps));       % set dP in FW compartment to 0
        %     else
        %         dp_FW = self.get_pressure_change_on_side(self.P_factor_FW,'FW');
        %     end
        %     if self.width_HW == 0
        %         dp_HW = zeros(length(self.y), length(self.time_steps));       % set dP in HW compartment to 0
        %     else
        %         dp_HW = self.get_pressure_change_on_side(self.P_factor_HW,'HW');
        %     end
        %     dp_fault = self.set_fault_pressure(dp_HW, dp_FW, self.dp_fault_mode);
        %     dp_fault = self.P_factor_fault' .* dp_fault;
        % end

        function p_fault = get_p_fault(self) 
            % get_p_fault Returns the pressure difference within the fault.
            % Output:
            %   p_fault - Pressure change within the fault
            % if reservoir compartment on either side has 0 width, set
            % pressures in that compartment to 0. 
            if self.width_FW == 0
                dp_FW = zeros(length(self.y), length(self.time_steps));       % set dP in FW compartment to 0
            else
                dp_FW = self.get_pressure_change_on_side(self.P_factor_FW,'FW');
            end
            if self.width_HW == 0
                dp_HW = zeros(length(self.y), length(self.time_steps));       % set dP in HW compartment to 0
            else
                dp_HW = self.get_pressure_change_on_side(self.P_factor_HW,'HW');
            end
            p0_FW = self.get_p0_on_side('FW');
            p0_HW = self.get_p0_on_side('HW');
            p_FW = p0_FW + dp_FW;
            p_HW = p0_HW + dp_HW;
            p0_fault = self.get_initial_pressure();
            % plot(p_FW(:,end),self.y, p_HW(:,end), self.y, p0_fault, self.y)
            % plot(p_FW(:,end),self.y, p_HW(:,end), self.y, p_fault, self.y)
            % dp_fault = self.set_fault_pressure(p_HW - p0_fault, p_FW - p0_fault, self.p_fault_mode);
            p_fault = self.set_fault_pressure(p_HW, p_FW, self.p_fault_mode);
            % set initial pressure again
            % otherwise different p0_fault_mode and p_fault_mode settings
            % lead to inconsistent state at t = 0;
            if self.time_steps(1) == 0
                p_fault(:,1) = p0_fault;
            end
            % p_fault = dp_fault + p0_fault
            % dp_fault = self.P_factor_fault' .* dp_fault;
        end


        function [fault_pressure] = set_fault_pressure(~, p_HW, p_FW, side_mode)
            % set_fault_pressure Sets the fault pressure based on HW and FW pressures.
            % Input:
            %   p_HW - Pressure in the hanging wall
            %   p_FW - Pressure in the footwall
            %   side_mode - Mode for setting the fault pressure
            % Output:
            %   fault_pressure - Calculated fault pressure
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
            elseif strcmp(side_mode, 'min_abs')
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

        function [p0_on_side] = get_p0_on_side(self, p_side)
            % get_p0_on_side returns pressure on HW or FW side
            % Input:
            % p_side - 'HW' or 'FW'
            % Output:
            % p0_on_side 

            if strcmp(p_side, 'FW')
                p0_on_side = self.get_FW_p0();
            elseif strcmp(p_side, 'HW')
                p0_on_side = self.get_HW_p0();
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end

        end

        function [p0_FW] = get_FW_p0(self)
            % get_FW_p0 return initial pressure in the footwall compartment
            % Output:
            %   p0_FW - Initial pressure in the footwall
            [next_to_FW, ~, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw);
            yy = self.y + self.depth_mid;    
            p0_FW = -(yy/1000).*self.p_grad + self.p_offset;       % [MPa] hydrostatic pressure
            top_FW_i = self.top_FW_i(self.y);                      % index where FW compartment starts (top)
            base_FW_i = self.base_FW_i(self.y);
            base_HW_i =  self.base_HW_i(self.y);
            % set overpressure w.r.t. hydrostatic gradient, in reservoir
            % and base
            p0_FW(top_FW_i:end) = p0_FW(top_FW_i:end) + self.p_over;  
            % set reservoir pressure gradient within the reservoir compartments
            if strcmp(self.p_res_mode, 'same')
                % if same, pressure is equal at the base of the deepest
                % compartment
                base_res_i = max(base_FW_i, base_HW_i );
                p0_FW(top_FW_i:base_res_i) = p0_FW(base_res_i) - (yy(top_FW_i:base_res_i) - yy(base_res_i))*self.p_grad_res/1000;
            else
                p0_FW(next_to_FW) = p0_FW(base_FW_i) - (yy(next_to_FW) - yy(base_FW_i))*self.p_grad_res/1000;
            end

        end

        function [p0_HW] = get_HW_p0(self)
            % get_HW_p0 return initial pressure in the footwall compartment
            % Output:
            %   p0_HW - Initial pressure in the footwall
            [~, next_to_HW, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw); 
            yy = self.y + self.depth_mid;    
            p0_HW = -(yy/1000).*self.p_grad + self.p_offset;       % [MPa] hydrostatic pressure
            top_HW_i = self.top_HW_i(self.y);                      % index where HW compartment starts (top)
            base_HW_i = self.base_HW_i(self.y);
            base_FW_i = self.base_FW_i(self.y);
            % set overpressure w.r.t. hydrostatic gradient, in reservoir
            % and base. Overpressure is defined at top reservoir
            % compartment
            p0_HW(top_HW_i:end) = p0_HW(top_HW_i:end) + self.p_over;  
            % set reservoir pressure gradient within the reservoir compartments
            if strcmp(self.p_res_mode, 'same')
                base_res_i = max(base_FW_i, base_HW_i );
                p0_HW(top_HW_i:base_res_i) = p0_HW(base_res_i) - (yy(top_HW_i:base_res_i) - yy(base_res_i))*self.p_grad_res/1000;
            else
                p0_HW(next_to_HW) = p0_HW(base_HW_i) - (yy(next_to_HW) - yy(base_HW_i))*self.p_grad_res/1000; 
            end
        end
       
        function  [dp_on_side] = get_pressure_change_on_side(self, p_factor, p_side)
            % get_pressure_change_on_side Calculates pressure change on a given side.
            % Input:
            %   p_factor - Pressure factor for the side
            %   p_side - Side indicator ('FW' or 'HW')
            % Output:
            %   dp_on_side - Pressure change on the specified side
            p_on_side = self.get_pressure_on_side(p_factor, p_side);
            p0_on_side = self.get_p0_on_side(p_side);
            dp_on_side = p_on_side - p0_on_side;
        end

        function  [dp_HW] = get_HW_pressure_change(self)
            % get_HW_pressure_change Returns pressure change in the hanging wall.
            % Output:
            %   dp_HW - Pressure change in the hanging wall
            p_HW = self.get_HW_pressure();
            p0_HW = self.get_HW_p0();
            dp_HW = p_HW - p0_HW;    
        end


        function  [dp_FW] = get_FW_pressure_change(self)
            % get_FW_pressure_change Returns pressure change in the footwall.
            % Output:
            %   dp_FW - Pressure change in the footwall
            p_FW = self.get_FW_pressure();
            p0_FW = self.get_FW_p0();
            dp_FW = p_FW - p0_FW;    
        end

        function [p_HW] = get_HW_pressure(self)
            % get_HW_pressure Returns pressure in the hanging wall.
            % Output:
            %   p_HW - Pressure in the hanging wall
            p_HW = self.get_pressure_on_side(self.P_factor_HW, 'HW');
        end

        function [p_FW] = get_FW_pressure(self)
            % get_FW_pressure Returns pressure in the footwall.
            % Output:
            %   p_FW - Pressure in the footwall
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
            p0_on_side = self.get_p0_on_side(p_side);

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
             
             p_on_side = dp_on_side + p0_on_side;
             if self.diffusion_P
                p_on_side = calc_dp_diffusion(self.y, y_top, y_base, self.time_steps, p_on_side, self.hyd_diffusivity);
             end
        end

       function a = get.p0(self)
            % get.p0 Returns the initial pressure.
            % Output:
            %   a - Initial pressure
            a =  self.get_initial_pressure();
       end
       
       function a = get.dp_fault(self)
            % get.dp_fault Returns the pressure difference across the fault.
            % Output:
            %   a - Pressure difference across the fault
            if contains(self.load_case,'P')
                a = self.get_p_fault() - self.get_initial_pressure;
                % a = self.get_dp_fault();
            else
                a = zeros(length(self.y), length(self.time_steps));
            end
       end

       function p = get.p(self)
            % get.p Returns the current pressure.
            % Output:
            %   p - Current pressure
            if contains(self.load_case,'P')
                p = self.get_p_fault();
            else
                p = repmat(self.p0, 1, length(self.time_steps));
            end
       end

        function self = reduce_steps(self, steps)
            % reduce_steps Reduces the number of steps in the simulation.
            % Input:
            %   steps - Number of steps to reduce to
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

        function [p_at_load_step] = get_p_at_load_step(self, load_step)
            % get_p_at_load_step Returns pressure at a specific load step.
            % Input:
            %   load_step - Load step to get pressure for
            % Output:
            %   p_at_load_step - Pressure at the specified load step
            % obtain the pressure at certain loadstep, where load_step is
            % between 1 and height of loadtable
            p_at_load_step = zeros(size(self.p, 1), 1);
            if load_step < 1 || load_step > size(self.p, 2)
                disp('Selected load step is outside calculation time range');
            else
                for i = 1 : size(self.p, 1)
                    x_ind = linspace(1, size(self.p, 2), size(self.p, 2));% indices of time, P, or T steps
                    p_at_load_step(i) = interp1(x_ind, self.p(i, :), load_step);
                end     
            end
        end

        function [dp_at_load_step] = get_dp_at_load_step(self, load_step)
            % get_dp_at_load_step Returns pressure difference at a specific load step.
            % Input:
            %   load_step - Load step to get pressure difference for
            % Output:
            %   dp_at_load_step - Pressure difference at the specified load step
            % obtain the pressure at certain loadstep, where load_step is
            % between 1 and height of loadtable
            dp_at_load_step = zeros(size(self.dp_fault, 1), 1);
            if load_step < 1 || load_step > size(self.dp, 2)
                disp('Selected load step is outside calculation time range');
            else
                for i = 1 : size(self.dp_fault, 1)
                    x_ind = linspace(1, size(self.dp_fault, 2), size(self.dp_fault, 2));    % indices of time, P, or T steps
                    dp_at_load_step(i) = interp1(x_ind, self.dp_fault(i, :), load_step);
                end     
            end
        end

        function check_if_property_exists(self, queried_property_name)
            % check_if_property_exists Checks if a property exists in the class.
            % Input:
            %   queried_property_name - Name of the property to check
            if ~ismember(queried_property_name, properties(self))
                error(['Queried property name is not a property of Class Pressure. ',newline,...
                    'Properties must be one of: ', strjoin(properties(self),', ')]);
            end

        end

        function [self] = update_properties(self, input)
            % update_properties Updates properties from input parameters.
            % Input:
            %   input - Input parameters (table or class)
            if istable(input)
                self = self.updatePropertiesFromTable(input);
            elseif isobject(input)
                self = self.updatePropertiesFromClass(input);
            end
        end

        function plot(self)
            % plot Plots the pressure and pressure change in the fault.
            subplot(1,2,1)
            plot(self.p(:,1));
            hold on
            plot(self.p(:,end));
            subplot(1,2,2)
            plot(self.dp_fault(:,end));
     
        end

        function self = updatePropertiesFromTable(self, inputTable)
            % updatePropertiesFromTable Updates properties from a table.
            % Input:
            %   inputTable - Table containing property values
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
            % updatePropertiesFromClass Updates properties from another class.
            % Input:
            %   inputClass - Class containing property values
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