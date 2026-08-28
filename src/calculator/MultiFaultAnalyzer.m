classdef MultiFaultAnalyzer < handle
    % MultiFaultAnalyzer handles multiple 2D fault cross-sections.
    %
    % The class stores one FaultAnalyzer object per fault and provides
    % convenience methods to:
    % - assign uniform or depth-dependent input parameters,
    % - set run settings across all faults,
    % - run all faults in serial or parallel,
    % - summarize and post-process outputs.
    %
    % Properties:
    %   faults - Array of FaultAnalyzer objects (one per fault)
    %   faultMetadata - Table with metadata per fault (ID, coordinates, etc.)
    %   faultSummary - Table summarizing run results per fault
    %   runDone - Logical flag indicating whether run() has completed
    %   parallel - Enables parallel execution when true
    %   printStatusOutput - Print per-fault progress output (default: true)
    %   printStatusEveryNFaults - Print every N faults (default: 1 = every fault)
    %
    % Dependent Properties:
    %   nFaults - Number of faults
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
        faults FaultAnalyzer      % array of FaultAnalyzer objects
        faultMetadata table         % table with custom meta data per fault (e.g. name, coordinates). ID is always included
        faultSummary table          % summary of fault results, e.g. reactivation & nucleation timestep, cff rate, slip length, etc. 
        runDone logical
        parallel = 1                % parallel calculation of faults stresses and slip
        printStatusOutput = true    % print per-fault progress during run
        printStatusEveryNFaults = 1 % print status every N faults (1 = every fault)
    end

    properties (Dependent)
        nFaults double
    end

    methods
        function self = MultiFaultAnalyzer()
            % MultiFaultCalculator Constructor to initialize the class with n_faults.
        end
        
        function self = initialize(self, nFaults, metadataTable)
        % Input:
            %   nFaults - Number of faults to be initialized
            % Optional:
            %   metadataTable   - table with additional metadata columns
            % construct the class with nFaults, assign ID in the
            % metadata table, and initialize the default FaultAnalyzer for each fault
            if nargin < 3 || isempty(metadataTable)
                metadataTable = table();
            end
            if ~istable(metadataTable)
                error('metadata_table must be a table');
            end

            self.faults = FaultAnalyzer.empty(0, 1);
            self.faultMetadata = table((1:nFaults)', 'VariableNames', {'ID'});
            if ~isempty(metadataTable)
                if height(metadataTable) ~= nFaults
                    error(['Metadata table height should match number of faults. #faults = ', num2str(nFaults), ...
                        ', #metadata rows = ', num2str(height(metadataTable))]);
                end
                self = self.addFaultMetadataAsTable(metadataTable);
            end
            for i = 1 : nFaults
                self.faults(i, 1) = FaultAnalyzer();
            end
            self.runDone = zeros(nFaults, 1);
        end

        function self = run(self)
            % run Runs the simulation for all faults.
            %
            % Each fault is handled in two steps that are both inside the
            % parallel loop so that Pressure, Temperature, and Green's
            % function construction (the expensive operations) are
            % distributed across workers:
            %
            %   extractInputs()                — builds plain-struct inputs
            %                                     (Pressure, Temperature, GF)
            %   computeStressAndNucleation() — static; stress + slip + nuc
            %
            % Results are then applied back serially via applyResults().
            % Using cell arrays for parfor input/output ensures MATLAB slices
            % exactly one FaultAnalyzer object per worker rather than
            % broadcasting the full typed array.
            all_faults  = self.faults;
            n           = self.nFaults;
            printStatus = self.printStatusOutput;
            printEveryN = max(1, round(self.printStatusEveryNFaults));

            % Pack faults into cell array for parfor slicing
            fault_cell = cell(n, 1);
            for i = 1 : n
                fault_cell{i} = all_faults(i);
            end

            result_cell = cell(n, 1);
            if self.parallel
                parfor i = 1 : n
                    inputs          = fault_cell{i}.extractInputs();
                    result_cell{i}  = FaultAnalyzer.computeStressAndNucleation(inputs);
                    if printStatus && (i == 1 || mod(i, printEveryN) == 0 || i == n)
                        fprintf('fault %d of %d\n', i, n);
                    end
                end
            else
                for i = 1 : n
                    inputs          = fault_cell{i}.extractInputs();
                    result_cell{i}  = FaultAnalyzer.computeStressAndNucleation(inputs);
                    if printStatus && (i == 1 || mod(i, printEveryN) == 0 || i == n)
                        fprintf('fault %d of %d\n', i, n);
                    end
                end
            end

            % Apply results serially
            for i = 1 : n
                all_faults(i) = all_faults(i).applyResults(result_cell{i});
            end
            self.faults       = all_faults;
            self.runDone      = true;
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
            nFaults = self.nFaults;
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

        function self = deactivateDepthDependentInputParameter(self, parameterName, indices)
            % deactivateDepthDependentInputParameter Switches selected
            % faults back to uniform-with-depth mode for one parameter.
            nFaults = self.nFaults;
            if nargin < 3 || isempty(indices)
                indices = 1:nFaults;
            end

            if ~isnumeric(indices) || any(indices < 1) || any(indices > nFaults) || any(mod(indices,1) ~= 0)
                error('indices must contain valid integer fault indices between 1 and n_faults');
            end
            indices = indices(:);

            for k = 1 : numel(indices)
                i = indices(k);
                self.faults(i) = self.faults(i).deactivateDepthDependentInputParameter(parameterName);
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
            if  (length(parameterValues) == self.nFaults) || isscalar(parameterValues)
                for i = 1 : self.nFaults
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
                if ~(height(inputTable) == self.nFaults)
                error(['Input table height should match number of faults on the fault',...
                    ' # of faults = ', num2str(self.nFaults), ' but # of table rows is ', num2str(height(inputTable))]);
            end
            
            % Get the list of properties of the fault input parameters 
            fault_input_props = properties(self.faults(1).faultParameterSpecs);
                        
            % Get the list of column names from the input table
            tableColumns = inputTable.Properties.VariableNames;
            
            % Loop through each property and update if there's a matching column in the table
            for j = 1:length(tableColumns)
                propName = tableColumns{j};
                if ismember(propName, fault_input_props)
                    for i = 1 : self.nFaults
                        self.faults(i) = self.faults(i).setInputParameter(propName, inputTable.(propName)(i));
                    end
                end
            end
        end

        function self = updateInputParameterFromMetadata(self)
            % updateInputParameterFromMetadata Updates input parameters from fault metadata.
            % Input values are applied per fault using the metadata columns.
            fault_input_props = properties(self.faults(1).faultParameterSpecs);
            tableColumns = self.faultMetadata.Properties.VariableNames;
            valid_props = intersect(tableColumns, fault_input_props, 'stable');
            nFaults = self.nFaults;
            for i = 1 : nFaults
                for j = 1 : length(valid_props)
                    propName = valid_props{j};
                    value = self.faultMetadata.(propName)(i);
                    if iscell(value)
                        value = value{1};
                    end
                    self.faults(i).setInputParameter(propName, value);
                end
                self.faults(i).realizationStale = true;
            end
        end

        function [inputParameters] = getInputParameter(self, parameterName)
            % getInputParameter Retrieves input parameters from all faults.
            % Input:
            %   parameter_name - Name of the parameter
            % Retrieve input parameters from all faults.
            % INPUT
            % parameter_name    string. parameter name, e.g. 'dip'
 
              inputParameters = zeros(self.nFaults, 1);
  %           if self.isValidInputParameterName(parameter_name)
                  for i = 1 : self.nFaults
                    inputParameters(i) = self.faults(i).getInputParameter(parameterName);
                end
   %         end
        end

        function [depths] = getDepth(self)
            % getDepth Returns the absolute depth vector for each fault.
            % Output:
            %   depths - Cell array with one absolute-depth vector per fault.
            depths = cell(self.nFaults, 1);
            for i = 1 : self.nFaults
                depths{i} = self.faults(i).getDepth();
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
                faultIndices = 1:self.nFaults;
            end

            if ~isnumeric(faultIndices) || any(mod(faultIndices, 1) ~= 0) || ...
                    any(faultIndices < 1) || any(faultIndices > self.nFaults)
                error('faultIndices must contain valid integer indices between 1 and nFaults');
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

        function [parameterValues] = getDepthDependentInputParameter(self, parameterName)
            % getDepthDependentInputParameter Retrieves depth-dependent
            % input parameter arrays from all faults.
            % Input:
            %   parameter_name - Name of the parameter
            % Output:
            %   parameterValues - cell array (n_faults x 1) with one
            %   depth-profile vector per fault

            n = self.nFaults;
            parameterValues = cell(n, 1);
            if n == 0
                return;
            end

            valid_input_parameter_names = properties(self.faults(1).faultParameterSpecs);
            if ~ismember(parameterName, valid_input_parameter_names)
                valid_names = [append(valid_input_parameter_names, repmat({', '}, length(valid_input_parameter_names), 1))];
                error(['input parameter name ', parameterName, ' not valid, should be one of ', valid_names{:}]);
            end

            for i = 1 : n
                parameterValues{i} = self.faults(i).getDepthDependentInputParameter(parameterName);
            end
        end

        function results = getResult(self, resultName, varargin)
            % getResult Return a result from one or more faults.
            %
            % For one fault, results is the value returned by FaultAnalyzer.
            % For multiple faults, results is a cell array in fault order.
            [faultIndices, allowStale] = self.parseResultQueryArguments(varargin{:});

            if isscalar(faultIndices)
                results = self.faults(faultIndices).getResult(resultName, 'AllowStale', allowStale);
                return;
            end

            results = cell(numel(faultIndices), 1);
            for i = 1:numel(faultIndices)
                results{i} = self.faults(faultIndices(i)).getResult(resultName, 'AllowStale', allowStale);
            end
        end

        function results = getResultAtLoadStep(self, resultName, loadStep, varargin)
            % getResultAtLoadStep Return a result at a load step for faults.
            [faultIndices, allowStale] = self.parseResultQueryArguments(varargin{:});
            results = self.collectFaultResults(@(fault) fault.getResultAtLoadStep(resultName, loadStep, 'AllowStale', allowStale), faultIndices);
        end

        function results = getResultAtY(self, resultName, yValue, varargin)
            % getResultAtY Return a result at a y value for faults.
            [faultIndices, allowStale] = self.parseResultQueryArguments(varargin{:});
            results = self.collectFaultResults(@(fault) fault.getResultAtY(resultName, yValue, 'AllowStale', allowStale), faultIndices);
        end

        function results = getResultAtDepth(self, resultName, depthValue, varargin)
            % getResultAtDepth Return a result at absolute depth for faults.
            [faultIndices, allowStale] = self.parseResultQueryArguments(varargin{:});
            results = self.collectFaultResults(@(fault) fault.getResultAtDepth(resultName, depthValue, 'AllowStale', allowStale), faultIndices);
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

            n = self.nFaults;
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
                self.faults(i).markResultsStale();
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
           if ~(size(loadTableArray,1) == 1 | size(loadTableArray,1) == self.nFaults)
               disp('ERROR: Input cell array of load tables must be n_faults x 1, or 1 x 1. Value not assigned');
           else
               for i = 1 : self.nFaults
                   loadTable = ensure_load_table_step_num(loadTableArray{min(i, size(loadTableArray, 1))});
                   self.faults(i).load_table = loadTable;
                   self.faults(i).markResultsStale();
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
            for i = 1 : self.nFaults
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
            if nargin < 2 || isempty(timeStepIndices)
                return;
            end

            if isscalar(timeStepIndices) && isnan(timeStepIndices)
                for i = 1 : self.nFaults
                    if ~self.faults(i).keepModelObjects
                        if isstruct(self.faults(i).faultResults)
                            if isfield(self.faults(i).faultResults, 'P'); self.faults(i).faultResults.P = []; end
                            if isfield(self.faults(i).faultResults, 'dP'); self.faults(i).faultResults.dP = []; end
                            if isfield(self.faults(i).faultResults, 'T'); self.faults(i).faultResults.T = []; end
                            if isfield(self.faults(i).faultResults, 'dT'); self.faults(i).faultResults.dT = []; end
                            if isfield(self.faults(i).faultResults, 'sne'); self.faults(i).faultResults.sne = []; end
                            if isfield(self.faults(i).faultResults, 'tau'); self.faults(i).faultResults.tau = []; end
                            if isfield(self.faults(i).faultResults, 'slip'); self.faults(i).faultResults.slip = []; end
                        end
                        self.faults(i).load_table = self.faults(i).load_table([],:);
                        continue;
                    end

                    if isobject(self.faults(i).stress{1})
                        self.faults(i).stress{1} = self.faults(i).stress{1}.reduce_steps(nan);
                    elseif isstruct(self.faults(i).stress{1})
                        if isfield(self.faults(i).stress{1}, 'sne'); self.faults(i).stress{1}.sne = []; end
                        if isfield(self.faults(i).stress{1}, 'tau'); self.faults(i).stress{1}.tau = []; end
                    end

                    if isobject(self.faults(i).temperature{1})
                        self.faults(i).temperature{1} = self.faults(i).temperature{1}.reduce_steps(nan);
                    elseif isstruct(self.faults(i).temperature{1})
                        if isfield(self.faults(i).temperature{1}, 'T'); self.faults(i).temperature{1}.T = []; end
                        if isfield(self.faults(i).temperature{1}, 'dT'); self.faults(i).temperature{1}.dT = []; end
                    end

                    if isobject(self.faults(i).pressure{1})
                        self.faults(i).pressure{1} = self.faults(i).pressure{1}.reduce_steps(nan);
                    elseif isstruct(self.faults(i).pressure{1})
                        if isfield(self.faults(i).pressure{1}, 'P'); self.faults(i).pressure{1}.P = []; end
                        if isfield(self.faults(i).pressure{1}, 'dP'); self.faults(i).pressure{1}.dP = []; end
                    end

                    if isobject(self.faults(i).slip{1})
                        self.faults(i).slip{1} = self.faults(i).slip{1}.reduce_steps(nan);
                    elseif isstruct(self.faults(i).slip{1})
                        if isfield(self.faults(i).slip{1}, 'slip'); self.faults(i).slip{1}.slip = []; end
                    end

                    if isstruct(self.faults(i).faultResults)
                        if isfield(self.faults(i).faultResults, 'P'); self.faults(i).faultResults.P = []; end
                        if isfield(self.faults(i).faultResults, 'dP'); self.faults(i).faultResults.dP = []; end
                        if isfield(self.faults(i).faultResults, 'T'); self.faults(i).faultResults.T = []; end
                        if isfield(self.faults(i).faultResults, 'dT'); self.faults(i).faultResults.dT = []; end
                        if isfield(self.faults(i).faultResults, 'sne'); self.faults(i).faultResults.sne = []; end
                        if isfield(self.faults(i).faultResults, 'tau'); self.faults(i).faultResults.tau = []; end
                        if isfield(self.faults(i).faultResults, 'slip'); self.faults(i).faultResults.slip = []; end
                    end
                    self.faults(i).load_table = self.faults(i).load_table([],:);
                end
                return;
            end

            if ~isnumeric(timeStepIndices) || any(~isfinite(timeStepIndices)) || any(mod(timeStepIndices,1) ~= 0)
                error('timeStepIndices must be a numeric vector of finite integers, or scalar NaN');
            end

            timeStepIndices = unique(timeStepIndices(:)');
            if min(timeStepIndices) < 1
                error('timeStepIndices must be >= 1');
            end

            for i = 1 : self.nFaults
                if ~self.faults(i).keepModelObjects
                    nSteps = size(self.faults(i).faultResults.sne, 2);
                else
                    nSteps = size(self.faults(i).stress{1}.sne, 2);
                end
                if max(timeStepIndices) > nSteps
                    error('timeStepIndices exceed available number of timesteps (%d) for fault %d', nSteps, i);
                end

                if ~self.faults(i).keepModelObjects
                    if isstruct(self.faults(i).faultResults)
                        if isfield(self.faults(i).faultResults, 'P'); self.faults(i).faultResults.P = self.faults(i).faultResults.P(:, timeStepIndices); end
                        if isfield(self.faults(i).faultResults, 'dP'); self.faults(i).faultResults.dP = self.faults(i).faultResults.dP(:, timeStepIndices); end
                        if isfield(self.faults(i).faultResults, 'T'); self.faults(i).faultResults.T = self.faults(i).faultResults.T(:, timeStepIndices); end
                        if isfield(self.faults(i).faultResults, 'dT'); self.faults(i).faultResults.dT = self.faults(i).faultResults.dT(:, timeStepIndices); end
                        if isfield(self.faults(i).faultResults, 'sne'); self.faults(i).faultResults.sne = self.faults(i).faultResults.sne(:, timeStepIndices); end
                        if isfield(self.faults(i).faultResults, 'tau'); self.faults(i).faultResults.tau = self.faults(i).faultResults.tau(:, timeStepIndices); end
                        if isfield(self.faults(i).faultResults, 'slip'); self.faults(i).faultResults.slip = self.faults(i).faultResults.slip(:, timeStepIndices); end
                    end
                    self.faults(i).load_table = self.faults(i).load_table(timeStepIndices,:);
                    continue;
                end

                if isobject(self.faults(i).stress{1})
                    self.faults(i).stress{1} = self.faults(i).stress{1}.reduce_steps(timeStepIndices);
                elseif isstruct(self.faults(i).stress{1})
                    if isfield(self.faults(i).stress{1}, 'sne'); self.faults(i).stress{1}.sne = self.faults(i).stress{1}.sne(:, timeStepIndices); end
                    if isfield(self.faults(i).stress{1}, 'tau'); self.faults(i).stress{1}.tau = self.faults(i).stress{1}.tau(:, timeStepIndices); end
                end

                if isobject(self.faults(i).temperature{1})
                    self.faults(i).temperature{1} = self.faults(i).temperature{1}.reduce_steps(timeStepIndices);
                elseif isstruct(self.faults(i).temperature{1})
                    if isfield(self.faults(i).temperature{1}, 'T'); self.faults(i).temperature{1}.T = self.faults(i).temperature{1}.T(:, timeStepIndices); end
                    if isfield(self.faults(i).temperature{1}, 'dT'); self.faults(i).temperature{1}.dT = self.faults(i).temperature{1}.dT(:, timeStepIndices); end
                end

                if isobject(self.faults(i).pressure{1})
                    self.faults(i).pressure{1} = self.faults(i).pressure{1}.reduce_steps(timeStepIndices);
                elseif isstruct(self.faults(i).pressure{1})
                    if isfield(self.faults(i).pressure{1}, 'P'); self.faults(i).pressure{1}.P = self.faults(i).pressure{1}.P(:, timeStepIndices); end
                    if isfield(self.faults(i).pressure{1}, 'dP'); self.faults(i).pressure{1}.dP = self.faults(i).pressure{1}.dP(:, timeStepIndices); end
                end

                if isobject(self.faults(i).slip{1})
                    self.faults(i).slip{1} = self.faults(i).slip{1}.reduce_steps(timeStepIndices);
                elseif isstruct(self.faults(i).slip{1})
                    if isfield(self.faults(i).slip{1}, 'slip'); self.faults(i).slip{1}.slip = self.faults(i).slip{1}.slip(:, timeStepIndices); end
                end

                if isstruct(self.faults(i).faultResults)
                    if isfield(self.faults(i).faultResults, 'P'); self.faults(i).faultResults.P = self.faults(i).faultResults.P(:, timeStepIndices); end
                    if isfield(self.faults(i).faultResults, 'dP'); self.faults(i).faultResults.dP = self.faults(i).faultResults.dP(:, timeStepIndices); end
                    if isfield(self.faults(i).faultResults, 'T'); self.faults(i).faultResults.T = self.faults(i).faultResults.T(:, timeStepIndices); end
                    if isfield(self.faults(i).faultResults, 'dT'); self.faults(i).faultResults.dT = self.faults(i).faultResults.dT(:, timeStepIndices); end
                    if isfield(self.faults(i).faultResults, 'sne'); self.faults(i).faultResults.sne = self.faults(i).faultResults.sne(:, timeStepIndices); end
                    if isfield(self.faults(i).faultResults, 'tau'); self.faults(i).faultResults.tau = self.faults(i).faultResults.tau(:, timeStepIndices); end
                    if isfield(self.faults(i).faultResults, 'slip'); self.faults(i).faultResults.slip = self.faults(i).faultResults.slip(:, timeStepIndices); end
                end

                self.faults(i).load_table = self.faults(i).load_table(timeStepIndices,:);
            end
        end


        function [depthMidValues] = getInputParameterAlongDepthMid(self, parameterName)
            % check whether the parsed input parameter name is valid
            [valid_input] = self.isValidInputParameterName(parameterName);
            depthMidValues = nan(1, self.nFaults);
            absolute_depths = self.getDepth();
            if valid_input
                for i = 1 : self.nFaults
                    parameter = self.faults(i).faultParameterSpecs.(parameterName);
                    if isnan(parameter.value_with_depth) | parameter.uniform_with_depth
                        depthMidValues(i) = parameter.value;
                    else
                        value_with_depth = parameter.value_with_depth;
                        depth = absolute_depths;
                        depth_mid = self.faults(i).faultParameterSpecs.depth_mid.value;
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
            reservoirBoundaries = array2table(zeros((self.nFaults),4),...
                'VariableNames', vars);
            for i = 1 : self.nFaults
                y = self.faults(i).y;
                for j = 1 : length(vars)
                    if isempty(self.faults(i).faultRealization)
                        self.faults(i).generateRealization();
                    end
                    reservoirBoundaries.(vars{j})(i) = self.faults(i).faultRealization.(['y_', vars{j}]) + self.faults(i).faultRealization.depth_mid;
                end
            end
        end

        function [minDepth, maxDepth] = getMinMaxDepth(self)
            % getMinMaxDepth Gets the shallowest and deepest point on the fault surface 
            absolute_depths = self.getDepth();
            for i = 1 : self.nFaults
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

        function [summary] = getResultsSummary(self)
            % getResultsSummary Gets the summary of results of individual
            % faults
            if ~isscalar(self.runDone) || ~self.runDone
                disp('Run not yet executed, empty summary');
                summary = table();
                return;
            end

            nFaults = self.nFaults;
            if nFaults == 0
                summary = table();
                return;
            end

            summaries = cell(nFaults, 1);
            for i = 1 : nFaults
                if ~istable(self.faults(i).faultSummary)
                    error('Fault %d does not contain a valid fault summary table.', i);
                end
                summaries{i} = self.faults(i).faultSummary;
            end
            summary = vertcat(summaries{:});
        end

        function self = regenerateFaultSummary(self, preserveCustomColumns)
            % regenerateFaultSummary Rebuilds faultSummary from per-fault
            % summaries and optionally preserves user-added columns.
            % Input:
            %   preserveCustomColumns - logical, default true
            if nargin < 2
                preserveCustomColumns = true;
            end

            refreshed = self.getResultsSummary();

            if ~preserveCustomColumns || isempty(refreshed)
                self.faultSummary = refreshed;
                return;
            end

            if isempty(self.faultSummary)
                self.faultSummary = refreshed;
                return;
            end

            if ~istable(self.faultSummary) || height(self.faultSummary) ~= height(refreshed)
                % Shapes differ; safest behavior is to use refreshed
                % summary only.
                self.faultSummary = refreshed;
                return;
            end

            merged = refreshed;
            refreshedNames = refreshed.Properties.VariableNames;
            existingNames = self.faultSummary.Properties.VariableNames;
            customNames = existingNames(~ismember(existingNames, refreshedNames));

            for i = 1:numel(customNames)
                cname = customNames{i};
                merged.(cname) = self.faultSummary.(cname);
            end

            self.faultSummary = merged;
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
            valid_field_names = fields(self.faults(1).faultParameterSpecs);
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

        function results = collectFaultResults(self, resultFunction, faultIndices)
            if isscalar(faultIndices)
                results = resultFunction(self.faults(faultIndices));
                return;
            end
            results = cell(numel(faultIndices), 1);
            for i = 1:numel(faultIndices)
                results{i} = resultFunction(self.faults(faultIndices(i)));
            end
        end

        function [faultIndices, allowStale] = parseResultQueryArguments(self, varargin)
            faultIndices = 1:self.nFaults;
            allowStale = false;
            optionArguments = varargin;
            if ~isempty(varargin) && isnumeric(varargin{1})
                faultIndices = varargin{1};
                optionArguments = varargin(2:end);
            end
            if ~isnumeric(faultIndices) || any(mod(faultIndices, 1) ~= 0) || ...
                    any(faultIndices < 1) || any(faultIndices > self.nFaults)
                error('faultIndices must contain valid integer indices between 1 and nFaults');
            end
            faultIndices = faultIndices(:);
            allowStale = false;
            if mod(numel(optionArguments), 2) ~= 0
                error('Options must be specified as name-value pairs.');
            end
            for i = 1:2:numel(optionArguments)
                if ~(ischar(optionArguments{i}) || (isstring(optionArguments{i}) && isscalar(optionArguments{i}))) || ...
                        ~strcmpi(char(optionArguments{i}), 'AllowStale')
                    error('Unknown option. Supported option: AllowStale.');
                end
                if ~(islogical(optionArguments{i + 1}) && isscalar(optionArguments{i + 1}))
                    error('AllowStale must be a logical scalar.');
                end
                allowStale = optionArguments{i + 1};
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
            if ismember(submittedName, valid_setting_names) & ~ismember(submittedName,{'faultParameterSpecs','input_parameters','load_table','y','realizationTable'})
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
            % TODO: still assumes timesteps are equal for all faults. 
            % Consider to constrain that all faults are ran with same nr of
            % load steps, and/or move this method to PantherAnalysis level
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

        function num = get.nFaults(self)
            % get.nFaults Gets the number of faults.
            % Output:
            %   num - Number of faults
            % number of faults in the object
            num = length(self.faults);
        end

    end
end





