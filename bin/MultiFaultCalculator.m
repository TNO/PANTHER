classdef MultiFaultCalculator
    % MultiFaultCalculator Class to perform calculations on multiple faults (pillars).
    % This class allows for different input and run settings for each pillar.
    %
    % Properties:
    %   pillars - Cell array of PANTHER input objects (1 ensemble per entry)
    %   pillar_info - Table with pillar custom meta_data (e.g. name, coordinates)
    %   pillar_results - Cell array to store results for each pillar
    %   result_summary - Table to summarize results
    %   run_done - Logical flag indicating if the run is completed
    %   parallel - Flag to enable parallel processing (default is 1)
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
        pillars cell        % cell array of PANTHER input objects (1 ensemble per entry, can be modified to multiple) 
        pillar_info table   % table with pillar custom meta_data (e.g. name, coordinates)
        pillar_results cell
        result_summary table
        run_done logical
        parallel = 1
    end

    properties (Dependent)
        n_pillars double
    end

    methods
        function self = MultiFaultCalculator(n_pillars)
            % MultiFaultCalculator Constructor to initialize the class with n_pillars.
            % Input:
            %   n_pillars - Number of pillars
            % construct the class with n_pillars, assign ID in the
            % information table
            self.pillars = cell(n_pillars, 1);
            self.pillar_info = table([1:n_pillars]','VariableNames',{'ID'});
            % initialize the default PANTHER input for each pillar
            for i = 1 : length(self.pillars)
                self.pillars{i} = PantherInput();
            end
        end

        function self = run(self)
            % run Runs the simulation for all fault pillars.
            all_pillars = self.pillars;
            n = self.n_pillars;
            self.pillar_results = cell(n, 1);
            if self.parallel
                parfor i = 1 : n
                    results{i,1} = panther(all_pillars{i});
                    disp([num2str(i),'/', num2str(n)]);
                end
            else
                results = cell(n,1);
                for i = 1 : n
                    results{i,1} = panther(all_pillars{i});
                    disp([num2str(i),'/', num2str(n)]);
                end
            end      
            self.pillar_results = results;
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
                self.pillar_info = [self.pillar_info, info_table_to_be_added];
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
                        disp(['depth dependent variable could not be set, size not equal to y, length is ', num2str(self.pillars{i}.y)]);
                    end
                end
            else
                disp(['variable ', parameter_name,' not assigned, check that length of input values equals number of pillars']);
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
            if self.is_valid_input_parameter_name(parameter_name) && (length(values) == length(self.pillars))
                for i = 1 : length(self.pillars)
                    self.pillars{i}.input_parameters.(parameter_name).(parameter_type) = values(i);
                end
            else
                disp(['variable ', parameter_name,' not assigned, check that length of input values equals number of pillars']);
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
            for i = 1 : length(self.pillar_results)
                nuc = new_nucleation_load_step;
                self.pillar_results{i}.stress{1} = self.pillar_results{i}.stress{1}.get_nucleation_stress(nuc);
            end
        end

        function self = reduce_output(self, time_step_indices)
            % reduce_output Reduces the output to given time step indices.
            % Input:
            %   time_step_indices - Indices of time steps to retain
            % return output only at give time step indices
            % provide nan if you don't want to store output (only reac and
            % nuc stresses are stored)
                for i = 1 : length(self.pillar_results)
                    if max(time_step_indices) < size(self.pillar_results{i}.stress{1}.sne, 2)  & ...
                        (min(time_step_indices) > 1)
                    self.pillar_results{i}.stress{1} = self.pillar_results{i}.stress{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.temperature{1} = self.pillar_results{i}.temperature{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.pressure{1} = self.pillar_results{i}.pressure{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.slip{1} = self.pillar_results{i}.slip{1}.reduce_steps(time_step_indices);
                        if ~isnan(time_step_indices)
                            self.pillar_results{i}.load_table = self.pillar_results{i}.load_table(time_step_indices,:);
                        end
                    end
                end
        end

        function self = add_info_from_closest_point(self, X, Y, Z_value, X_column, Y_column, new_column)
            % add_info_from_closest_point Adds meta data info based on nearest point.
            % Input:
            %   X - x-coordinates
            %   Y - y-coordinates
            %   Z_value - Spatial value
            %   X_column - Column name of the X_coordinate in pillar_info
            %   Y_column - Column name of the Y_coordinate in pillar_info
            %   new_column - Column name of Z_value
            % add meta data info based on nearest point 
            % INPUT
            % X             x-coordinates
            % Y             y-coordinates
            % Z_value       spatial value
            % X_column      string, column name of the X_coordinate in pillar_info
            % Y_column      string, column name of the X_coordinate in pillar_info
            % new_column    string, column name of Z_value
            i5 = floor(self.n_pillars/10);
            reverseStr = '';
            for i = 1 : length(self.pillars)
                xq = self.pillar_info.(X_column)(i);
                yq = self.pillar_info.(Y_column)(i);
                [closest_index] = self.nearest_rectangular_grid_coordinate(X,Y,xq, yq);
                if isempty(closest_index)
                    self.pillar_info.(new_column)(i) = Z_value(closest_index(1));
                    %self.display_progress(i, self.n_pillars, 20, ...
                    %    ['Progress assigning spatial data ', new_column,' : ']);
                    if rem(i, i5) == 0
                        percentDone = 100*i / self.n_pillars;
                        msg = sprintf('Progress assigning spatial data : %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                else
                    self.pillar_info.(new_column)(i) = Z_value(closest_index(1));
                end
            end
            sprintf(newline);
        end

        function [summary] = get_results_summary(self)
            % get_results_summary Gets the summary of results of individual
            % pillars
            if self.run_done
                summary = self.pillar_results{1}.summary;
                % concatenate summary tables of individual fault pillars
                for i = 1 : self.n_pillars - 1
                    summary = [summary; self.pillar_results{i}.summary];
                end
            else
                disp('Run not yet exectued, empty summary');
                summary = [];
            end
        end

        function [i_min] = nearest_rectangular_grid_coordinate(~, X_grid, Y_grid, X_query, Y_query)
            % nearest_rectangular_grid_coordinate Finds the nearest rectangular grid coordinate.
            % Input:
            %   X_grid - x-coordinates of the grid
            %   Y_grid - y-coordinates of the grid
            %   X_query - x-coordinate of the query point
            %   Y_query - y-coordinate of the query point
            % Output:
            %   i_min - Index of the nearest grid coordinate
            % only works for rectangular grids
            [~, i_min_x] = min(abs(X_query - X_grid));
            [~, i_min_y] = min(abs(Y_query - Y_grid));
            i_min = find((X_grid == X_grid(i_min_x)) & (Y_grid == Y_grid(i_min_y)));
            if isempty(i_min)
                distance = nan(size(X_grid));
                for i = 1 : length(X_grid)
                    distance(i) = pdist([X_query, Y_query; X_grid(i), Y_grid(i)]);
                    [~, i_min] = min(distance);
                end
                disp('The query point seems to be outside the grid, closest value was assinged');
            end
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
                disp(['Check length, or input parameter name not valid, should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        function [valid_name, value_type] = is_valid_setting_name(self, submitted_name)
            % is_valid_setting_name Checks if run setting name is valid.
            % Input:
            %   submitted_name - Name of the setting to validate
            % check if run setting name is valid
            valid_setting_names = fields(self.pillars{1});
            if ismember(submitted_name, valid_setting_names) & ~ismember(submitted_name,{'input_parameters','load_table','ensemble'})
                valid_name = true;
                if ismember(submitted_name,{'p_res_mode','p_fault_mode','dp_fault_mode','load_case','nucleation_criterion'})
                    value_type = 'char';
                else
                    value_type = 'double';
                end
            else
                valid_name = false;
                value_type = 'double';
                fields_cellstring = [append(valid_setting_names, repmat({', '},length(valid_setting_names),1))];
                disp(['Run setting name not valid, should be one of the following: ',...
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
    end
end
