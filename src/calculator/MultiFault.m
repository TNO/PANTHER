classdef MultiFault
    % MultiFault Class to perform calculations on multiple 2D cross-sections (faults)
    % along the same fault
    % This class allows for different input settings and run settings for each fault fault.
    % In addition, this class allows to specify additional metadata for the faults
    %
    % Properties:
    %   faults - Cell array of PANTHER input and results objects (1 ensemble per entry)
    %   fault_info - Table with custom meta_data (e.g. name, coordinates)
    %   fault_summary - Table to summarize results
    %   run_done - Logical flag indicating if the run is completed
    %   parallel - Flag to enable parallel processing (default is 1)
    %   suppress_faults_run_status_output - 1: do not display status update
    %   on number of faults that were processed
    %
    % Dependent Properties:
    %   n_faults - Number of faults
    %
    % Methods:
    %   MultiFault - Constructor to initialize the class with n_faults
    %   run - Runs the simulation for all fault faults
    %   add_fault_info_as_table - Adds meta data into the fault_info table
    %   set_depth_dependent_input_parameter - Sets depth-dependent input parameters
    %   set_input_parameter - Sets numeric input parameters
    %   update_input_parameter_from_table - Updates input parameters from a table
    %   get_input_value - Retrieves input parameters from the fault faults
    %   set_run_setting - Specifies run settings per fault
    %   get_minimum_nucleation_load_step - Gets the minimum nucleation load step
    %   overwrite_nucleation_stress - Overwrites nucleation stress
    %   reduce_output - Reduces the output to given time step indices
    %   add_info_from_closest_point - Adds meta data info based on nearest point
    %   get_results_summary - Gets the summary of results
    %   nearest_rectangular_grid_coordinate - Finds the nearest rectangular grid coordinate
    %   is_valid_input_parameter_name - Validates input parameter name
    %   is_valid_setting_name - Checks if run setting name is valid
    %   get.n_faults - Gets the number of faults

    properties
        faults PantherAnalysis  % cell array of PANTHER input & result objects (1 ensemble meber per entry, no stochastic analysis) 
        fault_metadata table        % table with custom meta_data per fault (e.g. name, coordinates)
        fault_summary table     % summary of fault results, e.g. reactivation & nucleation timestep, cff rate, slip length, etc. 
        run_done logical
        parallel = 1        % overrides parallel setting of individual faults
        suppress_fault_run_status_output = 0
    end

    properties (Dependent)
        n_faults double
    end

    methods
        function self = MultiFault()
            % MultiFaultCalculator Constructor to initialize the class with n_faults.
        end
        
        function self = initialize(self, n_faults)
        % Input:
            %   n_faults - Number of faults
            %   create_ensembles - optional logical flag to create initial ensembles on construction
            % construct the class with n_faults, assign ID in the
            % information table
            % initialize the default PANTHER model for each fault
            self.faults = PantherAnalysis.empty(0, 1);
            % Preallocate fault_info table
            self.fault_metadata = table((1:n_faults)', 'VariableNames', {'ID'});
            for i = 1 : n_faults
                % Construct PantherAnalysis without creating the ensemble (costly)
                self.faults(i, 1) = PantherAnalysis(false);
                % if create_ensembles
                %     % Generate ensembles only when explicitly requested
                %     self.faults(i, 1).generate_ensemble();
                % end
            end
            self.run_done = zeros(n_faults, 1);
        end

        function self = run(self)
            % run Runs the simulation for all fault faults.
            all_faults = self.faults;   % contains input objects for each fault
            n = self.n_faults;
            faults_updated_with_results = empty(n, 1);  % generate a separate output array to be able to use in parfor loop
            suppress_run_status_output = self.suppress_fault_run_status_output;
            if self.parallel
                parfor i = 1 : n
                    faults_updated_with_results(i,1) = panther(all_faults(i), self.suppress_fault_run_status_output);
                    if ~suppress_run_status_output
                        disp(['fault ', num2str(i),' of ', num2str(n)]);
                    end
                end
            else
                faults_updated_with_results = empty(n,1);
                for i = 1 : n
                    faults_updated_with_results(i,1) = panther(all_faults(i), self.suppress_fault_run_status_output);
                    if ~self.suppress_fault_run_status_output
                        disp(['fault ', num2str(i),' of ', num2str(n)]);
                    end
                end
            end      
            self.faults = faults_updated_with_results;        % update pantherinput objects with the results
            self.run_done = true;
            self.fault_summary = self.get_results_summary();  
        end

        function self = add_fault_metadata_as_table(self, info_table_to_be_added)
            % add_fault_info_as_table Adds meta data into the fault_info table.
            % Input:
            %   info_table_to_be_added - Table of height n_faults
            % adds meta data into the fault_info table
            if height(info_table_to_be_added) ~= height(self.fault_metadata)
                disp(['Cant append fault info, table size does not match. Height should be ',num2str(height(self.fault_metadata)) ]);
            else
                new_table_headers = info_table_to_be_added.Properties.VariableNames;
                [overlapping_headers, columns_in_fault_info] = ismember(new_table_headers, self.fault_metadata.Properties.VariableNames);
                if ~isempty(columns_in_fault_info)
                    self.fault_metadata(:, find(columns_in_fault_info))  = info_table_to_be_added(:, overlapping_headers);
                    self.fault_metadata = [self.fault_metadata, info_table_to_be_added(:, ~overlapping_headers)];
                else
                    self.fault_metadata = [self.fault_metadata, info_table_to_be_added];
                end
            end
        end

        function self = set_depth_dependent_input_parameter(self, parameter_name, parameter_values, indices)
            % set_depth_dependent_input_parameter Sets depth-dependent input parameters.
            % Input:
            %   parameter_name - Name of the parameter
            %   parameter_values - Cell array with one depth vector per selected fault,
            %                      or a single numeric vector when one fault index is selected.
            %   indices - Optional fault indices to update. If omitted,
            %             all faults are updated.
            % sets numeric input of depth-dependent Panther input parameters
            n_faults = length(self.faults);
            if nargin < 4 || isempty(indices)
                indices = 1:n_faults;
            end

            if ~isnumeric(indices) || any(indices < 1) || any(indices > n_faults) || any(mod(indices,1) ~= 0)
                error('indices must contain valid integer fault indices between 1 and n_faults');
            end
            indices = indices(:);

            % Convenience: for a single selected fault, allow passing the
            % depth profile directly as a numeric vector.
            if ~iscell(parameter_values)
                if numel(indices) == 1 && isnumeric(parameter_values) && isvector(parameter_values)
                    parameter_values = {parameter_values};
                else
                    error('Input depth dependent variable as a cell array (or numeric vector when assigning a single fault)');
                end
            end

            if isscalar(parameter_values)
                parameter_values = repmat(parameter_values, numel(indices), 1);
            elseif length(parameter_values) ~= numel(indices)
                error(['Depth dependent variable ', parameter_name, ' was not assigned, ', ...
                    'length of input cell array must equal number of selected faults']);
            end

            for k = 1 : numel(indices)
                i = indices(k);
                values_i = parameter_values{k};
                if length(values_i) == length(self.faults(i).y)
                    self.faults(i) = self.faults(i).set_depth_dependent_input_parameter(parameter_name, values_i);
                else
                    error(['Depth dependent variable ', parameter_name, ' could not be set for fault ', num2str(i), ...
                        ', size not equal to y, length of y is ', num2str(length(self.faults(i).y))]);
                end
            end
        end


        function self = set_input_parameter(self, parameter_name, parameter_values)
            % set_input_parameter Sets numeric input parameters.
            % Input:
            %   parameter_name - Name of the parameter
            %   parameter_values - Array of doubles length(faults), or
            %   single value
            % sets numeric input Panther input parameters
            % check the input parameter is the same length as n_faults
            if  (length(parameter_values) == self.n_faults) || isscalar(parameter_values)
                for i = 1 : length(self.faults)
                    if isscalar(parameter_values)
                        % same value assigned to all faults
                        self.faults(i).set_input_parameter(parameter_name, parameter_values);
                    else
                        self.faults(i).set_input_parameter(parameter_name, parameter_values(i));
                    end
                end
            else
                disp(['Variable ', parameter_name,' was not assigned, ',...
                ' length of input array does not equals number of faults']);
            end
        end
        

        function self = update_input_parameter_from_table(self, input_table)
            % update_input_parameter_from_table Updates input parameters from a table.
            % Input:
            %   input_table - Table containing property values
            %   parameter_type - Property type ('value', 'a', or 'b')

            % Check if height table matches number of fault
            if ~(height(input_table) == length(self.faults))
                error(['Input table height should match number of faults on the fault',...
                    ' # of faults = ', num2str(length(self.faults)), ' but # of table rows is ', num2str(height(input_table))]);
            end
            
            % Get the list of properties of the fault input parameters 
            fault_input_props = properties(self.faults(1).input_parameters);
                        
            % Get the list of column names from the input table
            tableColumns = input_table.Properties.VariableNames;
            
            % Loop through each property and update if there's a matching column in the table
            for j = 1:length(tableColumns)
                propName = tableColumns{j};
                if ismember(propName, fault_input_props)
                    for i = 1 : length(self.faults)
                        self.faults(i) = self.faults(i).set_input_parameter(propName, input_table.(propName)(i));
                    end
                end
            end
        end

        function self = update_input_parameter_from_metadata(self)
            % update_input_parameter_from_metadata Updates input parameters from fault metadata.
            % Input values are applied per fault using the metadata columns.
            fault_input_props = properties(self.faults(1).input_parameters);
            tableColumns = self.fault_info.Properties.VariableNames;
            valid_props = intersect(tableColumns, fault_input_props, 'stable');
            n_faults = self.n_faults;
            for i = 1 : n_faults
                for j = 1 : length(valid_props)
                    propName = valid_props{j};
                    value = self.fault_info.(propName)(i);
                    if iscell(value)
                        value = value{1};
                    end
                    self.faults(i).set_input_parameter(propName, value);
                end
                self.faults(i).ensemble_dirty = true;
            end
        end

        function [input_parameters] = get_input_parameter(self, parameter_name)
            % get_input_value Retrieves input parameters from the fault faults.
            % Input:
            %   parameter_name - Name of the parameter
            % retrieve input parameters from the fault faults
            % INPUT
            % parameter_name    string. parameter name, e.g. 'dip'
 
            input_parameters = zeros(self.n_faults, 1);
  %           if self.is_valid_input_parameter_name(parameter_name)
                 for i = 1 : self.n_faults
                    input_parameters(i) = self.faults(i).get_input_parameter(parameter_name);
                end
   %         end
        end

        function [parameter_values] = get_depth_dependent_input_parameter(self, parameter_name)
            % get_depth_dependent_input_parameter Retrieves depth-dependent
            % input parameter arrays from the fault faults.
            % Input:
            %   parameter_name - Name of the parameter
            % Output:
            %   parameter_values - cell array (n_faults x 1) with one
            %   depth-profile vector per fault

            n = self.n_faults;
            parameter_values = cell(n, 1);
            if n == 0
                return;
            end

            valid_input_parameter_names = properties(self.faults(1).input_parameters);
            if ~ismember(parameter_name, valid_input_parameter_names)
                valid_names = [append(valid_input_parameter_names, repmat({', '}, length(valid_input_parameter_names), 1))];
                error(['input parameter name ', parameter_name, ' not valid, should be one of ', valid_names{:}]);
            end

            for i = 1 : n
                parameter_values{i} = self.faults(i).get_depth_dependent_input_parameter(parameter_name);
            end
        end

        function self = set_run_setting(self, setting_name, setting_value)
            % set_run_setting Specifies run settings per fault.
            % Input:
            %   setting_name - Name of the setting to be applied
            %   setting_value - Cell array, array of floats, single cell, single float, or string
            % Specify run settings per fault.
            % INPUT
            % setting_name   - Name of the setting to be applied
            % setting_value  - Cell array, array of floats, single cell, single float, or string
        
            [is_valid, value_type] = self.is_valid_setting_name(setting_name);
            if ~is_valid
                error('Invalid setting name: %s', setting_name);
            end

            n = length(self.faults);
            % Fast path: if a scalar value is supplied, broadcast and assign
            % with minimal per-fault checks. For non-scalar inputs, convert
            % to a cell array of per-fault values and assign.
            if iscell(setting_value)
                if isscalar(setting_value)
                    assign_cells = repmat(setting_value(1), n, 1);
                elseif length(setting_value) == n
                    assign_cells = setting_value(:);
                else
                    error('Cell input must be scalar or length equal to number of faults');
                end
            elseif isnumeric(setting_value)
                if isscalar(setting_value)
                    assign_cells = repmat({setting_value}, n, 1);
                elseif length(setting_value) == n
                    assign_cells = num2cell(setting_value(:));
                else
                    error('Numeric input must be scalar or length equal to number of faults');
                end
            elseif ischar(setting_value) || isstring(setting_value)
                % broadcast string/char across faults
                assign_cells = repmat({char(setting_value)}, n, 1);
            else
                error('Unsupported setting_value type');
            end

            % Minimal assignment loop: no type checks inside the loop
            for i = 1 : n
                self.faults(i).(setting_name) = assign_cells{i};
            end
            
        end

        function self = set_load_tables(self, load_table_array)
            % load_table_array must be a cell array of size n_faults x 1,
            % or contain a single load table
            if ~iscell(load_table_array)
                disp('ERROR: Input cell array of load tables');
            end
           if ~(size(load_table_array,2) == 1)
               disp('ERROR: Input cell array of load tables must be n_faults x 1, or 1 x 1. Value not assigned');
           end
           if ~(size(load_table_array,1) == 1 | size(load_table_array,1) == length(self.faults))
               disp('ERROR: Input cell array of load tables must be n_faults x 1, or 1 x 1. Value not assigned');
           else
               for i = 1 : length(self.faults)
                   self.faults(i).load_table = load_table_array{i};
               end
           end
           
        end

        function [nuc_load_step] = get_minimum_nucleation_load_step(self)
            % get_minimum_nucleation_load_step Gets the minimum nucleation
            % load step over all the fault faults
            if self.run_done
                nuc_load_step = min(self.fault_summary.nucleation_load_step);
            else
                nuc_load_step = nan;
            end
        end

        function self = overwrite_nucleation_stress(self, new_nucleation_load_step)
            % overwrite_nucleation_stress Overwrites nucleation stress.
            % Input:
            %   new_nucleation_load_step - New nucleation load step
            for i = 1 : length(self.faults)
                nuc = new_nucleation_load_step;
                self.faults(i).stress{1} = self.faults(i).stress{1}.get_nucleation_stress(nuc);
            end
        end

        function self = reduce_output(self, time_step_indices)
            % reduce_output Reduces the output to given time step indices.
            % Input:
            % time_step_indices - Indices of time steps to retain
            % return output only at give time step indices
            % provide nan if you don't want to store output (only reac and
            % nuc stresses are stored)
                for i = 1 : length(self.faults)
                    if max(time_step_indices) < size(self.faults(i).stress{1}.sne, 2)  & ...
                        (min(time_step_indices) >= 1)
                    self.faults(i).stress{1} = self.faults(i).stress{1}.reduce_steps(time_step_indices);
                    self.faults(i).temperature{1} = self.faults(i).temperature{1}.reduce_steps(time_step_indices);
                    self.faults(i).pressure{1} = self.faults(i).pressure{1}.reduce_steps(time_step_indices);
                    self.faults(i).slip{1} = self.faults(i).slip{1}.reduce_steps(time_step_indices);
                        if ~isnan(time_step_indices)
                            self.faults(i).load_table = self.faults(i).load_table(time_step_indices,:);
                        end
                    end
                end
        end


        function [depth_mid_values] = get_input_parameter_along_depth_mid(self, parameter_name)
            % check whether the parsed input parameter name is valid
            [valid_input] = self.is_valid_input_parameter_name(parameter_name);
            depth_mid_values = nan(1, length(self.n_faults));
            absolute_depths = self.get_absolute_depths();
            if valid_input
                for i = 1 : length(self.faults)
                    parameter = self.faults(i).input_parameters.(parameter_name);
                    if isnan(parameter.value_with_depth) | parameter.uniform_with_depth
                        depth_mid_values(i) = parameter.value;
                    else
                        value_with_depth = parameter.value_with_depth;
                        depth = absolute_depths;
                        depth_mid = self.faults(i).input_parameters.depth_mid.value;
                        depth_mid_values(i) = interp1(depth, value_with_depth, depth_mid);    % should be the same as taking the middle element
                    end
                end
            else
                error([parameter_name, ' is not a valid input parameter name']);
            end
        end


        function [reservoir_boundaries] = get_top_base_reservoir(self)
            % get_top_base_reservoir  Gets the top and base of the footwall
            % and hanging wall reservoir compartments along the fault
            % strike
            % Ouput
            % reservoir_boundaries  (table) table of n_faultsx4, with
            % columns top_FW, base_FW, top_HW, base_HW
            vars = {'FW_top','FW_base','HW_top', 'HW_base'};
            reservoir_boundaries = array2table(zeros((self.n_faults),4),...
                'VariableNames', vars);
            for i = 1 : self.n_faults
                y = self.faults(i).y;
                for j = 1 : length(vars)
                    if isempty(self.faults(i).ensemble_members{1})
                        self.faults(i).generate_ensemble();
                    end
                    reservoir_boundaries.(vars{j})(i) = self.faults(i).ensemble_members{1}.(['y_', vars{j}]) + self.faults(i).ensemble_members{1}.depth_mid;
                end
            end
        end

        function [min_depth, max_depth] = get_min_max_depth(self)
            % get_min_max_depth Gets the shallowest and deepest point on the fault surface 
            absolute_depths = self.get_absolute_depths();
            for i = 1 : length(self.faults)
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
            % fault of the fault. 
            absolute_depths = cell(self.n_faults, 1);
            for i = 1 : length(self.faults)
                depth_mid = self.faults(i).input_parameters.depth_mid.value;
                absolute_depths{i} = self.faults(i).y + depth_mid;
            end
        end

        function [summary] = get_results_summary(self)
            % get_results_summary Gets the summary of results of individual
            % faults
            if self.run_done
                summary = self.faults(1).summary;
                % concatenate summary tables of individual fault faults
                for i = 1 : self.n_faults - 1
                    summary = [summary; self.faults(i).summary];
                end
            else
                disp('Run not yet exectued, empty summary');
                summary = [];
            end
        end

        function [valid_name] = is_valid_input_parameter_name(self, submitted_name, warning_on)
            % is_valid_input_parameter_name Validates input parameter name.
            % Input:
            %   submitted_name - Name of the parameter to validate
            % validate whether specified input parameter name is valid
            % Output:
            % valid_name: true or false
            if nargin < 3
                warning_on = true;
            end
            valid_field_names = fields(self.faults(1).input_parameters);
            if ismember(submitted_name, valid_field_names)
                valid_name = true;
            else
                valid_name = false;
                if warning_on
                    fields_cellstring = [append(valid_field_names, repmat({', '},length(valid_field_names),1))];
                    warning(['Submitted parameter name was ''', submitted_name,...
                        '''. Valid input parameter names are: ',...
                         [fields_cellstring{:}]]);
                end
            end
        end

        function [valid_name, output_category] = is_valid_output_name(self, submitted_name, warning_on)
            % is_valid_output_name Validates output parameter name.
            % Input:
            %   submitted_name - Name of the parameter to validate
            % validate whether specified output name is valid
            % Output:
            % valid_name: true or false
            if nargin < 3
                warning_on = true;
            end
            valid_field_names = {'P','dP','T','dT','sne','tau','slip'}';
            if ismember(submitted_name, valid_field_names)
                valid_name = true;
                if ismember(submitted_name, {'P','dP'})
                    output_category = 'pressure';
                elseif ismember(submitted_name, {'T','dT'})
                    output_category = 'temperature';
                elseif ismember(submitted_name, {'sne', 'tau'})
                    output_category = 'stress';
                elseif ismember(submitted_name, {'slip'})
                    output_category = 'slip';
                end
            else
                valid_name = false;
                if warning_on
                    fields_cellstring = [append(valid_field_names, repmat({', '},length(valid_field_names),1))];
                    warning(['Submitted parameter name was ''', submitted_name,...
                        '''. Valid output names are: ',...
                         [fields_cellstring{:}]]);
                end
                output_category = '';
            end
        end

        function [valid_name, value_type] = is_valid_setting_name(self, submitted_name)
            % is_valid_setting_name Checks if run setting name is valid.
            % Input:
            %   submitted_name - Name of the setting to validate
            % check if run setting name is valid
            valid_setting_names = fields(self.faults(1));
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
                disp(['Run setting name ''', submitted_name', ''' not valid, should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        function [valid_time_step] = is_valid_time_step(self, time_step)
            valid_time_step = false;
            if ~(time_step == floor(time_step))
                error(['Time step must be an integer between 1 and ',...
                    num2str(height(self.faults(1).load_table))]);
            elseif (time_step) > height(self.faults(1).load_table)
                error(['Specified time step ', num2str(time_step),...
                    ' exceeds number of time steps in the load table (', ...
                    num2str(height(self.faults(1).load_table)),')']);
            elseif (time_step) < 1
                error(['Time step must be an integer between 1 and ',...
                    num2str(height(self.faults(1).load_table))]);
            else 
                valid_time_step = true;
            end
        end

        function num = get.n_faults(self)
            % get.n_faults Gets the number of faults.
            % Output:
            %   num - Number of faults
            % number of faults in the object
            num = length(self.faults);
        end

    end
end
