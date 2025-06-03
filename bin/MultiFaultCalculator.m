classdef MultiFaultCalculator
    % MultiFaultCalculator Class to perform calculations on multiple 2D cross-sections (pillars)
    % along the same fault
    % This class allows for different input settings and run settings for each fault pillar.
    % In addition, this class allows to specify additional metadata for the pillars
    %
    % Properties:
    %   pillars - Cell array of PANTHER input and results objects (1 ensemble per entry)
    %   pillar_info - Table with pillar custom meta_data (e.g. name, coordinates)
    %   result_summary - Table to summarize results
    %   run_done - Logical flag indicating if the run is completed
    %   parallel - Flag to enable parallel processing (default is 1)
    %   suppress_pillar_run_status_output - 1: do not display status update
    %   on number of pillars that were processed
    %
    % Dependent Properties:
    %   n_pillars - Number of pillars
    %
    % Methods:
    %   MultiFaultCalculator - Constructor to initialize the class with n_pillars
    %   run - Runs the simulation for all fault pillars
    %   add_pillar_info_as_table - Adds meta data into the pillar_info table
    %   set_depth_dependent_input_parameter - Sets depth-dependent input parameters
    %   set_input_parameter - Sets numeric input parameters
    %   update_input_parameter_from_table - Updates input parameters from a table
    %   get_input_value - Retrieves input parameters from the fault pillars
    %   set_run_setting - Specifies run settings per pillar
    %   get_minimum_nucleation_load_step - Gets the minimum nucleation load step
    %   overwrite_nucleation_stress - Overwrites nucleation stress
    %   reduce_output - Reduces the output to given time step indices
    %   add_info_from_closest_point - Adds meta data info based on nearest point
    %   get_results_summary - Gets the summary of results
    %   nearest_rectangular_grid_coordinate - Finds the nearest rectangular grid coordinate
    %   is_valid_input_parameter_name - Validates input parameter name
    %   is_valid_setting_name - Checks if run setting name is valid
    %   get.n_pillars - Gets the number of pillars

    properties
        pillars cell        % cell array of PANTHER input & result objects (1 ensemble meber per entry, no stochastic analysis) 
        pillar_info table   % table with pillar custom meta_data (e.g. name, coordinates)
        result_summary table
        run_done logical
        parallel = 1        % overrides parallel setting of individual pillars
        suppress_pillar_run_status_output = 0
    end

    properties (Dependent)
        n_pillars double
        L double
    end

    methods
        function self = MultiFaultCalculator(n_pillars)
            % MultiFaultCalculator Constructor to initialize the class with n_pillars.
            % Input:
            %   n_pillars - Number of pillars
            % construct the class with n_pillars, assign ID in the
            % information table
            % initialize the default PANTHER model for each pillar
            self.pillars = cell(n_pillars, 1);
            for i = 1 : length(self.pillars)
                self.pillars{i} = PantherInput();
            end
            % add default information to the pillar_info table
            default_X_coordinates = linspace(0, n_pillars - 1, n_pillars)';     % horizontal coordinate X
            default_Y_coordinates = linspace(0, n_pillars - 1, n_pillars)';     % horizontal coordinate Y
            self.pillar_info = table((1:n_pillars)', default_X_coordinates, ...
                default_Y_coordinates,'VariableNames',{'ID','X','Y'});
            % suppresses output of every single pillar
            % instead, during the running of MultiFaultCalculator output
            % status will be given per fault
            self = self.set_run_setting('suppress_status_output', 1);
        end

        function self = run(self)
            % run Runs the simulation for all fault pillars.
            all_pillars = self.pillars;   % contains input objects for each pillar
            n = self.n_pillars;
            pillars_updated_with_results = cell(n, 1);  % generate a separate output array to be able to use in parfor loop
            suppress_run_status_output = self.suppress_pillar_run_status_output;
            if self.parallel
                parfor i = 1 : n
                    pillars_updated_with_results{i,1} = panther(all_pillars{i});
                    if ~suppress_run_status_output
                        disp(['Pillar ', num2str(i),' of ', num2str(n)]);
                    end
                end
            else
                pillars_updated_with_results = cell(n,1);
                for i = 1 : n
                    pillars_updated_with_results{i,1} = panther(all_pillars{i});
                    if ~self.suppress_pillar_run_status_output
                        disp(['Pillar ', num2str(i),' of ', num2str(n)]);
                    end
                end
            end      
            self.pillars = pillars_updated_with_results;        % update pantherinput objects with the results
            self.run_done = true;
            self.result_summary = self.get_results_summary();  
        end

        function self = add_pillar_info_as_table(self, info_table_to_be_added)
            % add_pillar_info_as_table Adds meta data into the pillar_info table.
            % Input:
            %   info_table_to_be_added - Table of height n_faults
            % adds meta data into the pillar_info table
            if height(info_table_to_be_added) ~= height(self.pillar_info)
                disp(['Cant append fault info, table size does not match. Height should be ',num2str(height(self.pillar_info)) ]);
            else
                new_table_headers = info_table_to_be_added.Properties.VariableNames;
                [overlapping_headers, columns_in_pillar_info] = ismember(new_table_headers, self.pillar_info.Properties.VariableNames);
                if ~isempty(columns_in_pillar_info)
                    self.pillar_info(:, find(columns_in_pillar_info))  = info_table_to_be_added(:, overlapping_headers);
                    self.pillar_info = [self.pillar_info, info_table_to_be_added(:, ~overlapping_headers)];
                else
                    self.pillar_info = [self.pillar_info, info_table_to_be_added];
                end
            end
        end

        function self = set_depth_dependent_input_parameter(self, parameter_name, values)
            % set_depth_dependent_input_parameter Sets depth-dependent input parameters.
            % Input:
            %   parameter_name - Name of the parameter
            %   values - Cell array of length n_pillars, containing arrays of doubles of length(y)
            % sets numeric input of depth-dependent Panther input parameters
            if ~iscell(values)
                error('Input depth dependent variable in a cell array of length n_pillars');
            end
            if self.is_valid_input_parameter_name(parameter_name) && (length(values) == length(self.pillars))
                for i = 1 : length(self.pillars)
                    if length(values{i}) == length(self.pillars{i}.y)
                        self.pillars{i}.input_parameters.(parameter_name).uniform_with_depth = 0;
                        self.pillars{i}.input_parameters.(parameter_name).value_with_depth = values{i};
                    else
                        disp(['Depth dependent variable could not be set, size not equal to y, length is ', num2str(self.pillars{i}.y)]);
                    end
                end
            else
                if ~self.is_valid_input_parameter_name(parameter_name)
                disp(['Depth dependent variable ', parameter_name,' was not assigned, ',...
                    ' wrong input name given']);
                elseif ~(length(values) == length(self.pillars))
                    disp(['Depth dependent variable ', parameter_name,' was not assigned, ',...
                    ' length of input array does not equals number of pillars']);
                end
            end
        end

        function self = set_input_parameter(self, parameter_name, values, parameter_type)
            % set_input_parameter Sets numeric input parameters.
            % Input:
            %   parameter_name - Name of the parameter
            %   values - Array of doubles
            %   parameter_type - Property type ('value', 'a', or 'b')
            % sets numeric input Panther input parameters
            if nargin < 4
                parameter_type = 'value';
            end
            valid_name = self.is_valid_input_parameter_name(parameter_name);
            if  valid_name && (length(values) == length(self.pillars))
                for i = 1 : length(self.pillars)
                    self.pillars{i}.input_parameters.(parameter_name).(parameter_type) = values(i);
                end
            else
                if ~valid_name
                    disp(['Variable ', parameter_name,' was not assigned, ',...
                    ' wrong input name given']);
                elseif ~(length(values) == length(self.pillars))
                    disp(['Variable ', parameter_name,' was not assigned, ',...
                    ' length of input array does not equals number of pillars']);
                end
            end
        end

        function self = update_input_parameter_from_table(self, input_table, parameter_type)
            % update_input_parameter_from_table Updates input parameters from a table.
            % Input:
            %   input_table - Table containing property values
            %   parameter_type - Property type ('value', 'a', or 'b')

            % Check if height table matches number of pillar
            if ~(height(input_table) == length(self.pillars))
                error(['Input table height should match number of pillars on the fault',...
                    ' # of pillars = ', num2str(length(self.pillars)), ' but # of table rows is ', num2str(height(input_table))]);
            end
 
            if nargin < 3
                parameter_type = 'value';
            end
            
            % Get the list of properties of the pillar input parameters 
            pillar_input_props = properties(self.pillars{1}.input_parameters);
                        
            % Get the list of column names from the input table
            tableColumns = input_table.Properties.VariableNames;
            
            % Loop through each property and update if there's a matching column in the table
            for j = 1:length(tableColumns)
                propName = tableColumns{j};
                if ismember(propName, pillar_input_props)
                    for i = 1 : length(self.pillars)
                        self.pillars{i}.input_parameters.(propName).(parameter_type) = input_table.(propName)(i);
                    end
                end
            end
        end

        function [input_values] = get_input_value(self, parameter_name, input_property)
            % get_input_value Retrieves input parameters from the fault pillars.
            % Input:
            %   parameter_name - Name of the parameter
            %   input_property - Property type ('value', 'a', or 'b')
            % retrieve input parameters from the fault pillars
            % INPUT
            % parameter_name    string. parameter name, e.g. 'dip'
            if nargin < 3
                input_property = 'value';
            end
            input_values = zeros(length(self.pillars), 1);
            if self.is_valid_input_parameter_name(parameter_name)
                 for i = 1 : length(self.pillars)
                    input_values(i) = self.pillars{i}.input_parameters.(parameter_name).(input_property);
                end
            end
        end

        function self = set_run_setting(self, setting_name, setting_value)
            % set_run_setting Specifies run settings per pillar.
            % Input:
            %   setting_name - Name of the setting to be applied
            %   setting_value - Cell array, array of floats, single cell, single float, or string
            % Specify run settings per pillar.
            % INPUT
            % setting_name   - Name of the setting to be applied
            % setting_value  - Cell array, array of floats, single cell, single float, or string
        
            [is_valid, value_type] = self.is_valid_setting_name(setting_name);
            
            if ~is_valid
                error('Invalid setting name: %s', setting_name);
            else
                valid_value = 1; % Initialize as valid
            end
            
            for i = 1 : length(self.pillars)
                if iscell(setting_value)
                    % Case 1: Cell array with the same length as pillars
                    if length(setting_value) == length(self.pillars)
                        assign_value = setting_value{i};
                    % Case 2: Single cell element (could be string, char array, or numeric)
                    elseif isscalar(setting_value)
                        assign_value = setting_value{1};
                    else
                        valid_value = 0; % Invalid case: cell array with incorrect length
                    end
                
                elseif isnumeric(setting_value)
                    % Case 3: Single numeric value
                    if isscalar(setting_value)
                        assign_value = setting_value(1);
                    % Case 4: Numeric array with the same length as pillars
                    elseif length(setting_value) == self.n_pillars
                        assign_value = setting_value(i);
                    else
                        valid_value = 0; % Invalid case: numeric array with incorrect length
                    end
                
                elseif ischar(setting_value) || isstring(setting_value)
                    % Case 5: String or character array
                    assign_value = setting_value;
                    
                else
                    valid_value = 0; % Invalid case: unsupported type
                end
                
                % Assign value if valid and type matches the expected type
                if valid_value && strcmp(value_type, class(assign_value))
                    self.pillars{i}.(setting_name) = assign_value;
                end
            end
            
            % Display a message if the value was not valid or the type didn't match
            if ~valid_value || ~strcmp(value_type, class(assign_value))
                disp('Specified setting type does not seem the right type or dimension, check');
            end
            
        end

        function self = set_load_tables(self, load_table_array)
            % load_table_array must be a cell array of size n_pillars x 1,
            % or contain a single load table
            if ~iscell(load_table_array)
                disp('ERROR: Input cell array of load tables');
            end
           if ~(size(load_table_array,2) == 1)
               disp('ERROR: Input cell array of load tables must be n_pillars x 1, or 1 x 1. Value not assigned');
           end
           if ~(size(load_table_array,1) == 1 | size(load_table_array,1) == length(self.pillars))
               disp('ERROR: Input cell array of load tables must be n_pillars x 1, or 1 x 1. Value not assigned');
           else
               for i = 1 : length(self.pillars)
                   self.pillars{i}.load_table = load_table_array{i};
               end
           end
           
        end

        function [nuc_load_step] = get_minimum_nucleation_load_step(self)
            % get_minimum_nucleation_load_step Gets the minimum nucleation
            % load step over all the fault pillars
            if self.run_done
                nuc_load_step = min(self.result_summary.nucleation_load_step);
            else
                nuc_load_step = nan;
            end
        end

        function self = overwrite_nucleation_stress(self, new_nucleation_load_step)
            % overwrite_nucleation_stress Overwrites nucleation stress.
            % Input:
            %   new_nucleation_load_step - New nucleation load step
            for i = 1 : length(self.pillars)
                nuc = new_nucleation_load_step;
                self.pillars{i}.stress{1} = self.pillars{i}.stress{1}.get_nucleation_stress(nuc);
            end
        end

        function self = reduce_output(self, time_step_indices)
            % reduce_output Reduces the output to given time step indices.
            % Input:
            %   time_step_indices - Indices of time steps to retain
            % return output only at give time step indices
            % provide nan if you don't want to store output (only reac and
            % nuc stresses are stored)
                for i = 1 : length(self.pillars)
                    if max(time_step_indices) < size(self.pillars{i}.stress{1}.sne, 2)  & ...
                        (min(time_step_indices) >= 1)
                    self.pillars{i}.stress{1} = self.pillars{i}.stress{1}.reduce_steps(time_step_indices);
                    self.pillars{i}.temperature{1} = self.pillars{i}.temperature{1}.reduce_steps(time_step_indices);
                    self.pillars{i}.pressure{1} = self.pillars{i}.pressure{1}.reduce_steps(time_step_indices);
                    self.pillars{i}.slip{1} = self.pillars{i}.slip{1}.reduce_steps(time_step_indices);
                        if ~isnan(time_step_indices)
                            self.pillars{i}.load_table = self.pillars{i}.load_table(time_step_indices,:);
                        end
                    end
                end
        end

        function self = add_info_from_closest_point(self, X_input, Y_input, Z_value, new_column_name, cutoff_distance, cutoff_value)
            % add_info_from_closest_point Adds meta data info based on nearest point.
            % Input:
            %   X_input - x-coordinates of data to append to the fault
            %   info. Can be a grid or a coordinate vector
            %   Y_input - y-coordinates of data to append to the fault info
            %   Can be a grid or a coordinate vector
            %   Z_value - Spatial value Can be a grid or a coordinate vector
            %   new_column_name - Column name of added Z_value in pillar_info
            %   cutoff distance - distance from fault pillar beyond which
            %   no value should be added
            %   cutoff_value - value to be specified for points too far
            %   away from the input coordinates
            if nargin < 7
                cutoff_value = NaN;
                if nargin < 6
                    cutoff_distance = 200;
                end
            end

            % Validation of input
            if ~isnumeric(cutoff_distance) || ~isnumeric(cutoff_value)
                error('Cutoff distance and value should be double');
                if max(size(cutoff_value)) > 1 || max(size(cutoff_distance)) > 1
                    error('Cutoff distance and value should be single numbers');
                end
            end
            if ~( (size(X_input, 1) == size(Y_input, 1)) & (size(Y_input, 1) == size(Z_value,1)) ) | ...
                ~( (size(X_input, 2) == size(Y_input, 2)) & (size(Y_input, 2) == size(Z_value,2)) )
                warning('Dimensions of X, Y input coordinates and input value are not equal, please check input');
                return
            end

            % parameters for status display
            i5 = floor(self.n_pillars/10);
            reverseStr = '';
            for i = 1 : length(self.pillars)
                if ~isvector(X_input)
                    X_input = X_input(:);
                    Y_input = Y_input(:);
                    Z_value = Z_value(:);
                end
                xq = self.pillar_info.X(i);
                yq = self.pillar_info.Y(i);
                [distance_to_query_point, closest_index] = self.index_of_nearest_coordinate(X_input,Y_input, xq, yq);
                if distance_to_query_point <= cutoff_distance
                    self.pillar_info.(new_column_name)(i) = Z_value(closest_index);
                else
                    self.pillar_info.(new_column_name)(i) = cutoff_value;
                end
                if rem(i, i5) == 0
                    percentDone = 100*i / self.n_pillars;
                    msg = sprintf('Progress assigning spatial data: %3.1f', percentDone);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
            end
            sprintf(newline);
        end

        function [L_grid, Z_grid, grid_value, grid_attributes] = get_fault_grid_for_input_parameter(self, parameter_name, depth_spacing)
            if nargin < 3
                depth_spacing = self.pillars{1}.dy;
            end
            % check whether the parsed input parameter name is valid
            self.is_valid_input_parameter_name(parameter_name);
            x_vector = self.L;
            [min_depth, max_depth] = self.get_min_max_depth();
            z_vector = (min_depth : depth_spacing : max_depth)';
            [L_grid, Z_grid] = meshgrid(x_vector, z_vector);
            absolute_depths = self.get_absolute_depths();
            grid_value = zeros(size(L_grid));
            for i = 1 : length(self.pillars)
                parameter = self.pillars{i}.input_parameters.(parameter_name);
                if isnan(parameter.value_with_depth) | parameter.uniform_with_depth
                    grid_value(:,i) = parameter.value;
                    grid_attributes.z_line(i,1) = parameter.value;
                    grid_attributes.uniform_value = true; 
                else
                    value_with_depth = parameter.value_with_depth;
                    depth = absolute_depths{i};
                    grid_value(:,i) = interp1(depth, value_with_depth, Z_grid(:,i));
                    grid_attributes.z_line(i,1) = nan;
                    grid_attributes.uniform_value = false; 
                end
            end
        end

        function [min_depth, max_depth] = get_min_max_depth(self)
            % get_min_max_depth Gets the shallowest and deepest point on the fault surface 
            absolute_depths = self.get_absolute_depths();
            for i = 1 : length(self.pillars)
                depth = absolute_depths{i};
                if i == 1
                    min_depth = min(depth);
                    max_depth = max(depth);
                else
                    min_depth = min(min_depth, min(depth));
                    max_depth = max(max_depth, max(depth));
                end
            end
        end

        function [absolute_depths] = get_absolute_depths(self)
            % get_absolute depths Gets the absolute depth range at each
            % pillar of the fault. 
            absolute_depths = cell(self.n_pillars, 1);
            for i = 1 : length(self.pillars)
                depth_mid = self.pillars{i}.input_parameters.depth_mid.value;
                absolute_depths{i} = self.pillars{i}.y + depth_mid;
            end
        end

        function [summary] = get_results_summary(self)
            % get_results_summary Gets the summary of results of individual
            % pillars
            if self.run_done
                summary = self.pillars{1}.summary;
                % concatenate summary tables of individual fault pillars
                for i = 1 : self.n_pillars - 1
                    summary = [summary; self.pillars{i}.summary];
                end
            else
                disp('Run not yet exectued, empty summary');
                summary = [];
            end
        end

        function [i_min, i_dist] = index_of_nearest_coordinate(~, X_vector, Y_vector, X_query, Y_query)
            % index_of_nearest_coordinate Finds the index of the nearest coordinate 
            % to x_query and y_query in the X and Y vectors.
            % Input:
            %   X_vector - x-coordinates
            %   Y_vector - y-coordinates 
            %   X_query - x-coordinate of the query point
            %   Y_query - y-coordinate of the query point
            % Output:
            %   i_min - Index of the nearest coordinate in X and Y vector
            %   i_dist - Distance of the nearest coordinate to query point
            distance_to_query_point = ((X_vector - X_query).^2 + (Y_vector - Y_query).^2).^0.5;
            [i_min, i_dist] = min(distance_to_query_point);
        end

        function [valid_name] = is_valid_input_parameter_name(self, submitted_name)
            % is_valid_input_parameter_name Validates input parameter name.
            % Input:
            %   submitted_name - Name of the parameter to validate
            % validate whether specified input parameter name is valid
            % Output:
            % valid_name: true or false
            valid_field_names = fields(self.pillars{1}.input_parameters);
            if ismember(submitted_name, valid_field_names)
                valid_name = true;
            else
                valid_name = false;
                fields_cellstring = [append(valid_field_names, repmat({', '},length(valid_field_names),1))];
                error(['Given input parameter name ', submitted_name,...
                    ' should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        function [valid_name, value_type] = is_valid_setting_name(self, submitted_name)
            % is_valid_setting_name Checks if run setting name is valid.
            % Input:
            %   submitted_name - Name of the setting to validate
            % check if run setting name is valid
            valid_setting_names = fields(self.pillars{1});
            if ismember(submitted_name, valid_setting_names) & ~ismember(submitted_name,{'input_parameters','load_table','y','ensemble'})
                valid_name = true;
                if ismember(submitted_name,{'P_res_mode','P_fault_mode','P0_fault_mode',...
                        'load_case','nucleation_criterion'})
                    value_type = 'char';
                else
                    value_type = 'double';
                end
            else
                valid_name = false;
                value_type = 'double';
                fields_cellstring = [append(valid_setting_names, repmat({', '},length(valid_setting_names),1))];
                disp(['Run setting name ', submitted_name', ' not valid, should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        function num = get.n_pillars(self)
            % get.n_pillars Gets the number of pillars.
            % Output:
            %   num - Number of pillars
            % number of faults in the object
            num = length(self.pillars);
        end

        function L = get.L(self)
            % get.L Gets the along-fault length based on the X and Y
            % coordinates of the pillars. Assumes the coordinates are in
            % the right order. 
            X = self.pillar_info.X;
            Y = self.pillar_info.Y;
            L = zeros(size(X));
            for i = 1 : length(X) - 1
                dx = X(i+1) - X(i);
                dy = Y(i+1) - Y(i);
                L(i+1) = (dx^2 + dy^2)^0.5;
            end
            L = cumsum(L);
        end

    end
end
