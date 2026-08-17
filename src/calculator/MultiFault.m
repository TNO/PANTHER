classdef MultiFault < handle
    % MultiFault handles multiple 2D fault cross-sections.
    %
    % The class stores one PantherAnalysis object per fault and provides
    % convenience methods to:
    % - assign uniform or depth-dependent input parameters,
    % - set run settings across all faults,
    % - run all faults in serial or parallel,
    % - summarize and post-process outputs.
    %
    % Properties:
    %   faults - Array of PantherAnalysis objects (one per fault)
    %   faultMetadata - Table with metadata per fault (ID, coordinates, etc.)
    %   faultSummary - Table summarizing run results per fault
    %   runDone - Logical flag indicating whether run() has completed
    %   parallel - Enables parallel execution when true
    %   suppress_fault_run_status_output - Suppress per-fault progress output
    %
    % Dependent Properties:
    %   n_faults - Number of faults
    %
    % Methods:
    %   initialize - Create default faults and optional metadata table
    %   run - Run all faults and populate faultSummary
    %   setInputParameter - Assign scalar input parameters across faults
    %   setDepthDependentInputParameter - Assign depth-dependent vectors
    %   getInputParameter - Retrieve scalar input parameters
    %   getDepthDependentInputParameter - Retrieve depth-dependent parameters
    %   setRunSetting - Set run settings for all faults
    %   getResultsSummary - Build summary table from fault results
    %

    properties
        faults PantherAnalysis      % array of PantherAnalysis objects
        faultMetadata table         % table with custom meta data per fault (e.g. name, coordinates). ID is always included
        faultSummary table          % summary of fault results, e.g. reactivation & nucleation timestep, cff rate, slip length, etc. 
        runDone logical
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
        
                 function self = initialize(self, nFaults, metadataTable)
        % Input:
            %   n_faults - Number of faults
            % Optional:
            %   metadata_table   - table with additional metadata columns
            % construct the class with n_faults, assign ID in the
            % metadata table, and initialize the default PANTHER model for each fault
            if nargin < 3 || isempty(metadataTable)
                metadataTable = table();
            end
            if ~istable(metadataTable)
                error('metadata_table must be a table');
            end

            self.faults = PantherAnalysis.empty(0, 1);
            self.faultMetadata = table((1:nFaults)', 'VariableNames', {'ID'});
            if ~isempty(metadataTable)
                if height(metadataTable) ~= nFaults
                    error(['Metadata table height should match number of faults. #faults = ', num2str(nFaults), ...
                        ', #metadata rows = ', num2str(height(metadataTable))]);
                end
                self = self.addFaultMetadataAsTable(metadataTable);
            end
            for i = 1 : nFaults
                self.faults(i, 1) = PantherAnalysis();
            end
            self.runDone = zeros(nFaults, 1);
        end

        function self = run(self)
            % run Runs the simulation for all faults.
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
            self.runDone = true;
            self.faultSummary = self.getResultsSummary();  
        end

        function self = addFaultMetadataAsTable(self, infoTableToBeAdded)
            % addFaultMetadataAsTable Adds metadata into the faultMetadata table.
            % Input:
            %   info_table_to_be_added - Table of height n_faults
            % adds meta data into the faultMetadata table
            if height(infoTableToBeAdded) ~= height(self.faultMetadata)
                disp(['Cant append fault info, table size does not match. Height should be ',num2str(height(self.faultMetadata)) ]);
            else
                new_table_headers = infoTableToBeAdded.Properties.VariableNames;
                [overlapping_headers, columns_in_fault_metadata] = ismember(new_table_headers, self.faultMetadata.Properties.VariableNames);
                if ~isempty(columns_in_fault_metadata)
                    self.faultMetadata(:, find(columns_in_fault_metadata))  = infoTableToBeAdded(:, overlapping_headers);
                    self.faultMetadata = [self.faultMetadata, infoTableToBeAdded(:, ~overlapping_headers)];
                else
                    self.faultMetadata = [self.faultMetadata, infoTableToBeAdded];
                end
            end
        end

        function self = setDepthDependentInputParameter(self, parameterName, parameterValues, indices)
            % setDepthDependentInputParameter Sets depth-dependent input parameters.
            % Input:
            %   parameter_name - Name of the parameter
            %   parameter_values - Cell array with one depth vector per selected fault,
            %                      or a single numeric vector when one fault index is selected.
            %   indices - Optional fault indices to update. If omitted,
            %             all faults are updated.
            % sets numeric input of depth-dependent Panther input parameters
            nFaults = self.n_faults;
            if nargin < 4 || isempty(indices)
                indices = 1:nFaults;
            end

            if ~isnumeric(indices) || any(indices < 1) || any(indices > nFaults) || any(mod(indices,1) ~= 0)
                error('indices must contain valid integer fault indices between 1 and n_faults');
            end
            indices = indices(:);

            % Convenience: for a single selected fault, allow passing the
            % depth profile directly as a numeric vector.
            if ~iscell(parameterValues)
                if isscalar(indices) && isnumeric(parameterValues) && isvector(parameterValues)
                    parameterValues = {parameterValues};
                else
                    error('Input depth dependent variable as a cell array (or numeric vector when assigning a single fault)');
                end
            end

            if isscalar(parameterValues)
                parameterValues = repmat(parameterValues, numel(indices), 1);
            elseif length(parameterValues) ~= numel(indices)
                error(['Depth dependent variable ', parameterName, ' was not assigned, ', ...
                    'length of input cell array must equal number of selected faults']);
            end

            for k = 1 : numel(indices)
                i = indices(k);
                values_i = parameterValues{k};
                if length(values_i) == length(self.faults(i).y)
                    self.faults(i) = self.faults(i).setDepthDependentInputParameter(parameterName, values_i);
                else
                    error(['Depth dependent variable ', parameterName, ' could not be set for fault ', num2str(i), ...
                        ', size not equal to y, length of y is ', num2str(length(self.faults(i).y))]);
                end
            end
        end


        function self = setInputParameter(self, parameterName, parameterValues)
            % setInputParameter Sets numeric input parameters.
            % Input:
            %   parameter_name - Name of the parameter
            %   parameter_values - Array of doubles length(faults), or
            %   single value
            % sets numeric input Panther input parameters
            % check the input parameter is the same length as n_faults
            if  (length(parameterValues) == self.n_faults) || isscalar(parameterValues)
                for i = 1 : self.n_faults
                    if isscalar(parameterValues)
                        % same value assigned to all faults
                        self.faults(i).setInputParameter(parameterName, parameterValues);
                    else
                        self.faults(i).setInputParameter(parameterName, parameterValues(i));
                    end
                end
            else
                disp(['Variable ', parameterName,' was not assigned, ',...
                ' length of input array does not equals number of faults']);
            end
        end
        

        function self = updateInputParameterFromTable(self, inputTable)
            % updateInputParameterFromTable Updates input parameters from a table.
            % Input:
            %   input_table - Table containing property values
            %   parameter_type - Property type ('value', 'a', or 'b')

            % Check if height table matches number of fault
                if ~(height(inputTable) == self.n_faults)
                error(['Input table height should match number of faults on the fault',...
                    ' # of faults = ', num2str(self.n_faults), ' but # of table rows is ', num2str(height(inputTable))]);
            end
            
            % Get the list of properties of the fault input parameters 
            fault_input_props = properties(self.faults(1).input_parameters);
                        
            % Get the list of column names from the input table
            tableColumns = inputTable.Properties.VariableNames;
            
            % Loop through each property and update if there's a matching column in the table
            for j = 1:length(tableColumns)
                propName = tableColumns{j};
                if ismember(propName, fault_input_props)
                    for i = 1 : self.n_faults
                        self.faults(i) = self.faults(i).setInputParameter(propName, inputTable.(propName)(i));
                    end
                end
            end
        end

        function self = updateInputParameterFromMetadata(self)
            % updateInputParameterFromMetadata Updates input parameters from fault metadata.
            % Input values are applied per fault using the metadata columns.
            fault_input_props = properties(self.faults(1).input_parameters);
            tableColumns = self.faultMetadata.Properties.VariableNames;
            valid_props = intersect(tableColumns, fault_input_props, 'stable');
            n_faults = self.n_faults;
            for i = 1 : n_faults
                for j = 1 : length(valid_props)
                    propName = valid_props{j};
                    value = self.faultMetadata.(propName)(i);
                    if iscell(value)
                        value = value{1};
                    end
                    self.faults(i).setInputParameter(propName, value);
                end
                self.faults(i).ensemble_dirty = true;
            end
        end

        function [inputParameters] = getInputParameter(self, parameterName)
            % getInputParameter Retrieves input parameters from all faults.
            % Input:
            %   parameter_name - Name of the parameter
            % Retrieve input parameters from all faults.
            % INPUT
            % parameter_name    string. parameter name, e.g. 'dip'
 
            inputParameters = zeros(self.n_faults, 1);
  %           if self.isValidInputParameterName(parameter_name)
                 for i = 1 : self.n_faults
                    inputParameters(i) = self.faults(i).getInputParameter(parameterName);
                end
   %         end
        end

        function [parameterValues] = getDepthDependentInputParameter(self, parameterName)
            % getDepthDependentInputParameter Retrieves depth-dependent
            % input parameter arrays from all faults.
            % Input:
            %   parameter_name - Name of the parameter
            % Output:
            %   parameterValues - cell array (n_faults x 1) with one
            %   depth-profile vector per fault

            n = self.n_faults;
            parameterValues = cell(n, 1);
            if n == 0
                return;
            end

            valid_input_parameter_names = properties(self.faults(1).input_parameters);
            if ~ismember(parameterName, valid_input_parameter_names)
                valid_names = [append(valid_input_parameter_names, repmat({', '}, length(valid_input_parameter_names), 1))];
                error(['input parameter name ', parameterName, ' not valid, should be one of ', valid_names{:}]);
            end

            for i = 1 : n
                parameterValues{i} = self.faults(i).getDepthDependentInputParameter(parameterName);
            end
        end

        function self = setRunSetting(self, settingName, settingValue)
            % setRunSetting Specifies run settings per fault.
            % Input:
            %   setting_name - Name of the setting to be applied
            %   setting_value - Cell array, array of floats, single cell, single float, or string
            % Specify run settings per fault.
            % INPUT
            % setting_name   - Name of the setting to be applied
            % setting_value  - Cell array, array of floats, single cell, single float, or string
        
            [is_valid, value_type] = self.isValidSettingName(settingName);
            if ~is_valid
                error('Invalid setting name: %s', settingName);
            end

            n = self.n_faults;
            % Fast path: if a scalar value is supplied, broadcast and assign
            % with minimal per-fault checks. For non-scalar inputs, convert
            % to a cell array of per-fault values and assign.
            if iscell(settingValue)
                if isscalar(settingValue)
                    assign_cells = repmat(settingValue(1), n, 1);
                elseif length(settingValue) == n
                    assign_cells = settingValue(:);
                else
                    error('Cell input must be scalar or length equal to number of faults');
                end
            elseif isnumeric(settingValue)
                if isscalar(settingValue)
                    assign_cells = repmat({settingValue}, n, 1);
                elseif length(settingValue) == n
                    assign_cells = num2cell(settingValue(:));
                else
                    error('Numeric input must be scalar or length equal to number of faults');
                end
            elseif ischar(settingValue) || isstring(settingValue)
                % broadcast string/char across faults
                assign_cells = repmat({char(settingValue)}, n, 1);
            else
                error('Unsupported setting_value type');
            end

            % Minimal assignment loop: no type checks inside the loop
            for i = 1 : n
                self.faults(i).(settingName) = assign_cells{i};
            end
            
        end

        function self = setLoadTables(self, loadTableArray)
            % load_table_array must be a cell array of size n_faults x 1,
            % or contain a single load table
            if ~iscell(loadTableArray)
                disp('ERROR: Input cell array of load tables');
            end
           if ~(size(loadTableArray,2) == 1)
               disp('ERROR: Input cell array of load tables must be n_faults x 1, or 1 x 1. Value not assigned');
           end
           if ~(size(loadTableArray,1) == 1 | size(loadTableArray,1) == self.n_faults)
               disp('ERROR: Input cell array of load tables must be n_faults x 1, or 1 x 1. Value not assigned');
           else
               for i = 1 : self.n_faults
                   self.faults(i).load_table = loadTableArray{i};
               end
           end
           
        end

        function [nucLoadStep] = getMinimumNucleationLoadStep(self)
            % getMinimumNucleationLoadStep Gets the minimum nucleation
            % load step over all faults.
            if self.runDone
                nucLoadStep = min(self.faultSummary.nucleation_load_step);
            else
                nucLoadStep = nan;
            end
        end

        function self = overwriteNucleationStress(self, newNucleationLoadStep)
            % overwriteNucleationStress Overwrites nucleation stress.
            % Input:
            %   new_nucleation_load_step - New nucleation load step
            for i = 1 : self.n_faults
                nuc = newNucleationLoadStep;
                self.faults(i).stress{1} = self.faults(i).stress{1}.get_nucleation_stress(nuc);
            end
        end

        function self = reduceOutput(self, timeStepIndices)
            % reduceOutput Reduces the output to given time step indices.
            % Input:
            % time_step_indices - Indices of time steps to retain
            % return output only at give time step indices
            % provide nan if you don't want to store output (only reac and
            % nuc stresses are stored)
                for i = 1 : self.n_faults
                    if max(timeStepIndices) < size(self.faults(i).stress{1}.sne, 2)  & ...
                        (min(timeStepIndices) >= 1)
                    self.faults(i).stress{1} = self.faults(i).stress{1}.reduce_steps(timeStepIndices);
                    self.faults(i).temperature{1} = self.faults(i).temperature{1}.reduce_steps(timeStepIndices);
                    self.faults(i).pressure{1} = self.faults(i).pressure{1}.reduce_steps(timeStepIndices);
                    self.faults(i).slip{1} = self.faults(i).slip{1}.reduce_steps(timeStepIndices);
                        if ~isnan(timeStepIndices)
                            self.faults(i).load_table = self.faults(i).load_table(timeStepIndices,:);
                        end
                    end
                end
        end


        function [depthMidValues] = getInputParameterAlongDepthMid(self, parameterName)
            % check whether the parsed input parameter name is valid
            [valid_input] = self.isValidInputParameterName(parameterName);
            depthMidValues = nan(1, length(self.n_faults));
            absolute_depths = self.getAbsoluteDepths();
            if valid_input
                for i = 1 : self.n_faults
                    parameter = self.faults(i).input_parameters.(parameterName);
                    if isnan(parameter.value_with_depth) | parameter.uniform_with_depth
                        depthMidValues(i) = parameter.value;
                    else
                        value_with_depth = parameter.value_with_depth;
                        depth = absolute_depths;
                        depth_mid = self.faults(i).input_parameters.depth_mid.value;
                        depthMidValues(i) = interp1(depth, value_with_depth, depth_mid);    % should be the same as taking the middle element
                    end
                end
            else
                error([parameterName, ' is not a valid input parameter name']);
            end
        end


        function [reservoirBoundaries] = getTopBaseReservoir(self)
            % getTopBaseReservoir  Gets the top and base of the footwall
            % and hanging wall reservoir compartments along the fault
            % strike
            % Ouput
            % reservoirBoundaries  (table) table of n_faultsx4, with
            % columns top_FW, base_FW, top_HW, base_HW
            vars = {'FW_top','FW_base','HW_top', 'HW_base'};
            reservoirBoundaries = array2table(zeros((self.n_faults),4),...
                'VariableNames', vars);
            for i = 1 : self.n_faults
                y = self.faults(i).y;
                for j = 1 : length(vars)
                    if isempty(self.faults(i).ensemble_members{1})
                        self.faults(i).generate_ensemble();
                    end
                    reservoirBoundaries.(vars{j})(i) = self.faults(i).ensemble_members{1}.(['y_', vars{j}]) + self.faults(i).ensemble_members{1}.depth_mid;
                end
            end
        end

        function [minDepth, maxDepth] = getMinMaxDepth(self)
            % getMinMaxDepth Gets the shallowest and deepest point on the fault surface 
            absolute_depths = self.getAbsoluteDepths();
            for i = 1 : self.n_faults
                depth = absolute_depths{i};
                if i == 1
                    minDepth = min(depth);
                    maxDepth = max(depth);
                else
                    minDepth = min(minDepth, min(depth));
                    maxDepth = max(maxDepth, max(depth));
                end
            end
        end

        function [absoluteDepths] = getAbsoluteDepths(self)
            % getAbsoluteDepths Gets the absolute depth range for each
            % fault.
            absoluteDepths = cell(self.n_faults, 1);
            for i = 1 : self.n_faults
                depth_mid = self.faults(i).input_parameters.depth_mid.value;
                absoluteDepths{i} = self.faults(i).y + depth_mid;
            end
        end

        function [LDip, dLDip] = getLAlongDip(self, faultIndices)
            % getLAlongDip Wrapper around PantherAnalysis.getLAlongDip for one or more faults.
            % Input:
            %   faultIndices - optional fault indices (default: all faults)
            % Output:
            %   LDip  - along-fault length vector (single fault) or cell array
            %   dLDip - spacing scalar/vector (single fault) or cell array
            if nargin < 2 || isempty(faultIndices)
                faultIndices = 1:self.n_faults;
            end

            if ~isnumeric(faultIndices) || any(mod(faultIndices, 1) ~= 0) || ...
                    any(faultIndices < 1) || any(faultIndices > self.n_faults)
                error('faultIndices must contain valid integer indices between 1 and n_faults');
            end
            faultIndices = faultIndices(:);

            if isscalar(faultIndices)
                [LDip, dLDip] = self.faults(faultIndices).getLAlongDip();
                return;
            end

            LDip = cell(numel(faultIndices), 1);
            dLDip = cell(numel(faultIndices), 1);
            for k = 1 : numel(faultIndices)
                [LDip{k}, dLDip{k}] = self.faults(faultIndices(k)).getLAlongDip();
            end
        end

        function [summary] = getResultsSummary(self)
            % getResultsSummary Gets the summary of results of individual
            % faults
            if self.runDone
                summary = self.faults(1).summary;
                % Concatenate summary tables from individual faults.
                for i = 1 : self.n_faults - 1
                    summary = [summary; self.faults(i).summary];
                end
            else
                disp('Run not yet exectued, empty summary');
                summary = [];
            end
        end

        function [validName] = isValidInputParameterName(self, submittedName, warningOn)
            % isValidInputParameterName Validates input parameter name.
            % Input:
            %   submitted_name - Name of the parameter to validate
            % validate whether specified input parameter name is valid
            % Output:
            % valid_name: true or false
            if nargin < 3
                warningOn = true;
            end
            valid_field_names = fields(self.faults(1).input_parameters);
            if ismember(submittedName, valid_field_names)
                validName = true;
            else
                validName = false;
                if warningOn
                    fields_cellstring = [append(valid_field_names, repmat({', '},length(valid_field_names),1))];
                    warning(['Submitted parameter name was ''', submittedName,...
                        '''. Valid input parameter names are: ',...
                         [fields_cellstring{:}]]);
                end
            end
        end

        function [validName, outputCategory] = isValidOutputName(self, submittedName, warningOn)
            % isValidOutputName Validates output parameter name.
            % Input:
            %   submitted_name - Name of the parameter to validate
            % validate whether specified output name is valid
            % Output:
            % valid_name: true or false
            if nargin < 3
                warningOn = true;
            end
            valid_field_names = {'P','dP','T','dT','sne','tau','slip'}';
            if ismember(submittedName, valid_field_names)
                validName = true;
                if ismember(submittedName, {'P','dP'})
                    outputCategory = 'pressure';
                elseif ismember(submittedName, {'T','dT'})
                    outputCategory = 'temperature';
                elseif ismember(submittedName, {'sne', 'tau'})
                    outputCategory = 'stress';
                elseif ismember(submittedName, {'slip'})
                    outputCategory = 'slip';
                end
            else
                validName = false;
                if warningOn
                    fields_cellstring = [append(valid_field_names, repmat({', '},length(valid_field_names),1))];
                    warning(['Submitted parameter name was ''', submittedName,...
                        '''. Valid output names are: ',...
                         [fields_cellstring{:}]]);
                end
                outputCategory = '';
            end
        end

        function [validName, valueType] = isValidSettingName(self, submittedName)
            % isValidSettingName Checks if run setting name is valid.
            % Input:
            %   submitted_name - Name of the setting to validate
            % check if run setting name is valid
            valid_setting_names = fields(self.faults(1));
            if ismember(submittedName, valid_setting_names) & ~ismember(submittedName,{'input_parameters','load_table','y','ensemble'})
                validName = true;
                if ismember(submittedName,{'P_res_mode','P_fault_mode','P0_fault_mode',...
                        'load_case','nucleation_criterion'})
                    valueType = 'char';
                else
                    valueType = 'double';
                end
            else
                validName = false;
                valueType = 'double';
                fields_cellstring = [append(valid_setting_names, repmat({', '},length(valid_setting_names),1))];
                disp(['Run setting name ''', submittedName', ''' not valid, should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        function [validTimeStep] = isValidTimeStep(self, timeStep)
            validTimeStep = false;
            if ~(timeStep == floor(timeStep))
                error(['Time step must be an integer between 1 and ',...
                    num2str(height(self.faults(1).load_table))]);
            elseif (timeStep) > height(self.faults(1).load_table)
                error(['Specified time step ', num2str(timeStep),...
                    ' exceeds number of time steps in the load table (', ...
                    num2str(height(self.faults(1).load_table)),')']);
            elseif (timeStep) < 1
                error(['Time step must be an integer between 1 and ',...
                    num2str(height(self.faults(1).load_table))]);
            else 
                validTimeStep = true;
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
