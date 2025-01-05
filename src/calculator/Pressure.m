classdef (HandleCompatible) Pressure < ModelGeometry & FaultMesh
    % PantherPressure Class to represent the pressure calculations in the geological model.
    % A simplified pressure profile is prescribed, with uniform pressure
    % change in the reservoir and difusion to the seal and base (optional)
    % This class extends ModelGeometry and FaultMesh to include properties and methods
    % for calculating pressure changes in the model, including initial
    % pressure, and pressure (changes) during load steps in the fault, as
    % well as methods to retrieve pressures in the HW and FW compartments
    %
    % Properties:
    %   hyd_diffusivity - Hydraulic diffusivity
    %   time_steps - Time steps for the simulation
    %   P_steps - Pressure change steps
    %   P_factor_HW - Pressure factor for the hanging wall
    %   P_factor_FW - Pressure factor for the footwall
    %   P_factor_fault - Pressure factor for the fault
    %   P0_fault_mode - Mode to set initial fault pressure
    %   P_fault_mode - Mode to set fault pressure during loading
    %   P_res_mode - Mode to set initial reservoir pressures
    %   diffusion_P - Logical flag for diffusion
    %   P_grad - Pressure gradient
    %   P_grad_res - Reservoir pressure gradient
    %   P_offset - Pressure offset
    %   P_over - Overpressure
    %   load_case - Load case identifier
    %
    % Dependent Properties:
    %   P0 - [MPa] Initial pressure in the fault
    %   P - [MPa] Pressure in the fault during load steps 
    %   dP - [MPa] Pressure change in the fault during load steps
    %
    % Methods:
    %   PantherPressure - Constructor to initialize the pressure model
    %   get_initial_pressure - Returns the initial pressure in the model
    %   get_dP_fault - Returns the pressure difference across the fault
    %   set_fault_pressure - Sets the fault pressure based on HW and FW pressures
    %   get_initial_pressure_on_side - Returns the HW or FW initial pressure 
    %   get_P0_HW - Returns initial pressure in the hanging wall
    %   get_P0_FW - Returns initial pressure in the footwall
    %   get_dP_on_side - Calculates pressure change on a given side
    %   get_dP_HW - Returns pressure change in the hanging wall
    %   get_dP_FW - Returns pressure change in the footwall
    %   get_P_HW - Returns pressure in the hanging wall
    %   get_P_FW - Returns pressure in the footwall
    %   get_P_on_side - Calculates pressure on a given side
    %   reduce_steps - Reduces the number of steps in the simulation
    %   get_P_at_load_step - Returns fault pressure at a specific load step
    %   get_dP_at_load_step - Returns fault pressure difference at a specific load step
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
        P0_fault_mode {mustBeMember(P0_fault_mode,{'max','min','mean','FW','HW'})} = 'max';
        P_fault_mode {mustBeMember(P_fault_mode,{'max','max_abs','min', 'min_abs','mean','FW','HW'})} = 'min'; 
        P_res_mode {mustBeMember(P_res_mode, {'same','different'})} = 'same'
        diffusion_P logical = false
        P_grad (:,1) double {mustBePositive} = 10.2
        P_grad_res (:,1) double {mustBePositive} = 10.2
        P_offset (:,1) double  = 0
        P_over (:,1) double  = 0
        load_case = 'P'
    end

    properties (Access=private)
    end

    properties (Dependent)
        P0          % Initial pressure in the fault
        P           % Pressure in the fault during load steps
        dP          % Pressure change in the fault during load steps
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

       
        function [P0] = get_initial_pressure(self)
            % get_initial_pressure Returns the initial fault pressure
            % Output:
            %   P0 - Initial pressure in the fault
            P0_FW = self.get_P0_FW();
            P0_HW = self.get_P0_HW();
            P0 = self.set_fault_pressure(P0_HW, P0_FW, self.P0_fault_mode);
        end


        function P_fault = get_P_fault(self) 
            % get_p_fault Returns the pressure difference within the fault.
            % Output:
            %   P_fault - Pressure change within the fault
            % if reservoir compartment on either side has 0 width, set
            % pressures in that compartment to 0. 
            if self.width_FW == 0
                dP_FW = zeros(length(self.y), length(self.time_steps));       % set dP in FW compartment to 0
            else
                dP_FW = self.get_dP_on_side(self.P_factor_FW,'FW');
            end
            if self.width_HW == 0
                dP_HW = zeros(length(self.y), length(self.time_steps));       % set dP in HW compartment to 0
            else
                dP_HW = self.get_dP_on_side(self.P_factor_HW,'HW');
            end
            P0_FW = self.get_P0_on_side('FW');
            P0_HW = self.get_P0_on_side('HW');
            P_FW = P0_FW + dP_FW;
            P_HW = P0_HW + dP_HW;
            P0_fault = self.get_initial_pressure();
            P_fault = self.set_fault_pressure(P_HW, P_FW, self.P_fault_mode);
            % set initial pressure again to ensure consistent state at t = 0;
            if self.time_steps(1) == 0
                P_fault(:,1) = P0_fault;
            end
        end


        function [fault_pressure] = set_fault_pressure(~, P_HW, P_FW, side_mode)
            % set_fault_pressure Sets the fault pressure based on HW and FW pressures.
            % Input:
            %   p_HW - Pressure in the hanging wall
            %   p_FW - Pressure in the footwall
            %   side_mode - Mode for setting the fault pressure
            % Output:
            %   fault_pressure - Calculated fault pressure
            fault_pressure = zeros(size(P_HW));
            if strcmp(side_mode, 'max')
                fault_pressure = max(P_FW, P_HW);
            elseif strcmp(side_mode, 'max_abs')
                for i = 1 : size(P_FW,2)
                    [~, ind] = max([abs(P_FW(:,i)), abs(P_HW(:,i))],[],2);
                    fault_pressure(ind == 1, i) = P_FW(ind == 1, i);
                    fault_pressure(ind == 2, i) = P_HW(ind == 2, i);
                end
            elseif strcmp(side_mode, 'min')
                fault_pressure = min(P_FW, P_HW);
            elseif strcmp(side_mode, 'min_abs')
                for i = 1 : size(P_FW,2)
                    [~, ind] = min([abs(P_FW(:,i)), abs(P_HW(:,i))],[],2);
                    fault_pressure(ind == 1, i) = P_FW(ind == 1, i);
                    fault_pressure(ind == 2, i) = P_HW(ind == 2, i);
                end
            elseif strcmp(side_mode, 'mean')
                for i = 1 : size(P_FW,2)
                    fault_pressure(:,i) = mean([P_FW(:,i), P_HW(:,i)],2);
                end
            elseif strcmp(side_mode, 'FW')
                fault_pressure = P_FW;
            elseif strcmp(side_mode, 'HW')
                fault_pressure = P_HW;
            end
        end

        function [P0_on_side] = get_P0_on_side(self, P_side)
            % get_p0_on_side returns pressure on HW or FW side
            % Input:
            % p_side - 'HW' or 'FW'
            % Output:
            % p0_on_side 

            if strcmp(P_side, 'FW')
                P0_on_side = self.get_P0_FW();
            elseif strcmp(P_side, 'HW')
                P0_on_side = self.get_P0_HW();
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end

        end

        function [P0_FW] = get_P0_FW(self)
            % get_P0_FW return initial pressure in the footwall compartment
            % Output:
            %   P0_FW - Initial pressure in the footwall
            [next_to_FW, ~, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw);
            yy = self.y + self.depth_mid;    
            P0_FW = -(yy/1000).*self.P_grad + self.P_offset;       % [MPa] hydrostatic pressure
            top_FW_i = self.top_FW_i(self.y);                      % index where FW compartment starts (top)
            base_FW_i = self.base_FW_i(self.y);
            base_HW_i =  self.base_HW_i(self.y);
            % set overpressure w.r.t. hydrostatic gradient, in reservoir
            % and base
            P0_FW(top_FW_i:end) = P0_FW(top_FW_i:end) + self.P_over;  
            % set reservoir pressure gradient within the reservoir compartments
            if strcmp(self.P_res_mode, 'same')
                % if same, pressure is equal at the base of the deepest
                % compartment
                base_res_i = max(base_FW_i, base_HW_i );
                P0_FW(top_FW_i:base_res_i) = P0_FW(base_res_i) - (yy(top_FW_i:base_res_i) - yy(base_res_i))*self.P_grad_res/1000;
            else
                P0_FW(next_to_FW) = P0_FW(base_FW_i) - (yy(next_to_FW) - yy(base_FW_i))*self.P_grad_res/1000;
            end

        end

        function [P0_HW] = get_P0_HW(self)
            % get_P0_HW return initial pressure in the footwall compartment
            % Output:
            %   P0_HW - Initial pressure in the footwall
            [~, next_to_HW, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw); 
            yy = self.y + self.depth_mid;    
            P0_HW = -(yy/1000).*self.P_grad + self.P_offset;       % [MPa] hydrostatic pressure
            top_HW_i = self.top_HW_i(self.y);                      % index where HW compartment starts (top)
            base_HW_i = self.base_HW_i(self.y);
            base_FW_i = self.base_FW_i(self.y);
            % set overpressure w.r.t. hydrostatic gradient, in reservoir
            % and base. Overpressure is defined at top reservoir
            % compartment
            P0_HW(top_HW_i:end) = P0_HW(top_HW_i:end) + self.P_over;  
            % set reservoir pressure gradient within the reservoir compartments
            if strcmp(self.P_res_mode, 'same')
                base_res_i = max(base_FW_i, base_HW_i );
                P0_HW(top_HW_i:base_res_i) = P0_HW(base_res_i) - (yy(top_HW_i:base_res_i) - yy(base_res_i))*self.P_grad_res/1000;
            else
                P0_HW(next_to_HW) = P0_HW(base_HW_i) - (yy(next_to_HW) - yy(base_HW_i))*self.P_grad_res/1000; 
            end
        end
       
        function  [dP_on_side] = get_dP_on_side(self, P_factor, P_side)
            % get_dP_on_side Calculates pressure change on a given side.
            % Input:
            %   P_factor - Pressure factor for the side
            %   P_side - Side indicator ('FW' or 'HW')
            % Output:
            %   dP_on_side - Pressure change on the specified side
            P_on_side = self.get_P_on_side(P_factor, P_side);
            P0_on_side = self.get_P0_on_side(P_side);
            dP_on_side = P_on_side - P0_on_side;
        end

        function  [dP_HW] = get_dP_HW(self)
            % get_dP_HW Returns pressure change in the hanging wall.
            % Output:
            %   dP_HW - Pressure change in the hanging wall
            P_HW = self.get_P_HW();
            P0_HW = self.get_P0_HW();
            dP_HW = P_HW - P0_HW;    
        end


        function  [dP_FW] = get_dP_FW(self)
            % get_dP_FW Returns pressure change in the footwall.
            % Output:
            %   dP_FW - Pressure change in the footwall
            P_FW = self.get_P_FW();
            P0_FW = self.get_P0_FW();
            dP_FW = P_FW - P0_FW;    
        end

        function [P_HW] = get_P_HW(self)
            % get_P_HW Returns pressure in the hanging wall.
            % Output:
            %   P_HW - Pressure in the hanging wall
            P_HW = self.get_P_on_side(self.P_factor_HW, 'HW');
        end

        function [P_FW] = get_P_FW(self)
            % get_P_FW Returns pressure in the footwall.
            % Output:
            %   P_FW - Pressure in the footwall
            P_FW = self.get_P_on_side(self.P_factor_FW, 'FW');
        end


        function [P_on_side] = get_P_on_side(self, P_factor, P_side)
            % get_P_on_side Method to calculate pressure on a given side (FW or HW).
            %
            % This method calculates the pressure in a given reservoir compartment (FW or HW) based on
            % the input pressure factor for that compartment and side indicator. It also accounts for diffusion
            % if specified.
            %
            % Parameters:
            %   P_factor - Pressure factor for the side.
            %   P_side - Indicator for reservoir compartment ('FW' or 'HW').
            %
            % Returns:
            %   P_on_side - Calculated pressure on the specified side.

            % Validate side indicator
            if strcmp(P_side, 'FW')
                [next_to_p_side, ~, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw);
            elseif strcmp(P_side, 'HW')
                [~, next_to_p_side, ~] = is_adjacent_to_reservoir(self.y, self.thick, self.throw);
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end
            P0_on_side = self.get_P0_on_side(P_side);

            % Get top and base y-coordinates for the reservoir compartment
            y_top = self.get_top_y(P_side);
            y_base = self.get_base_y(P_side);

            % Ensure dp_unit is a column vector
            dP_unit = double(next_to_p_side);       
            if ~iscolumn(dP_unit)
                 dP_unit = dP_unit';
            end

            % Ensure P_steps and p_factor are row vectors
            P_factor = P_factor(:)';
             
             dP_on_side = dP_unit * self.P_steps' .* P_factor;
             
             P_on_side = dP_on_side + P0_on_side;
             if self.diffusion_P
                P_on_side = calc_dp_diffusion(self.y, y_top, y_base, self.time_steps, P_on_side, self.hyd_diffusivity);
             end
        end

       function P0 = get.P0(self)
            % get.P0 Returns the initial pressure.
            % Output:
            %   P0 - Initial pressure
            P0 =  self.get_initial_pressure();
       end
       
       function dP = get.dP(self)
            % get.dP_fault Returns the pressure difference across the fault.
            % Output:
            %   dP - Pressure difference across the fault
            if contains(self.load_case,'P')
                dP = self.get_P_fault() - self.get_initial_pressure;
            else
                dP = zeros(length(self.y), length(self.time_steps));
            end
       end

       function P = get.P(self)
            % get.P Returns the pressure in the fault.
            % Output:
            %   P - Current pressure
            if contains(self.load_case,'P')
                P = self.get_P_fault();
            else
                P = repmat(self.P0, 1, length(self.time_steps));
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

        function [P_at_load_step] = get_P_at_load_step(self, load_step)
            % get_p_at_load_step Returns pressure at a specific load step.
            % Input:
            %   load_step - Load step to get pressure for
            % Output:
            %   P_at_load_step - Pressure at the specified load step
            % obtain the pressure at certain loadstep, where load_step is
            % between 1 and height of loadtable
            P_at_load_step = zeros(size(self.P, 1), 1);
            if load_step < 1 || load_step > size(self.P, 2)
                disp('Selected load step is outside calculation time range');
            else
                for i = 1 : size(self.P, 1)
                    x_ind = linspace(1, size(self.P, 2), size(self.P, 2));% indices of time, P, or T steps
                    P_at_load_step(i) = interp1(x_ind, self.P(i, :), load_step);
                end     
            end
        end

        function [dP_at_load_step] = get_dP_at_load_step(self, load_step)
            % get_dp_at_load_step Returns pressure difference at a specific load step.
            % Input:
            %   load_step - Load step to get pressure difference for
            % Output:
            %   dp_at_load_step - Pressure difference at the specified load step
            % obtain the pressure at certain loadstep, where load_step is
            % between 1 and height of loadtable
            dP_at_load_step = zeros(size(self.dP, 1), 1);
            if load_step < 1 || load_step > size(self.dp, 2)
                disp('Selected load step is outside calculation time range');
            else
                for i = 1 : size(self.dP, 1)
                    x_ind = linspace(1, size(self.dP, 2), size(self.dP, 2));    % indices of time, P, or T steps
                    dP_at_load_step(i) = interp1(x_ind, self.dP(i, :), load_step);
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
            plot(self.P(:,1));
            hold on
            plot(self.P(:,end));
            subplot(1,2,2)
            plot(self.dP(:,end));
     
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