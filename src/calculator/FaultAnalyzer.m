classdef (HandleCompatible) FaultAnalyzer < FaultMesh 
    % Object that initializes input, sets run and save settings for Panther, and
    % contains the results

    properties
        faultParameterSpecs                         % {FaultParameterList object} containing fault parameter specifications
        faultRealization = []                       % {FaultRealization object} for this fault; empty until generateRealization() is called
        faultResults struct                         % [struct] lightweight plain-array results (sne, tau, P, slip, etc.); always populated after run
        faultSummary table                          % [table] per-run scalar summary (reactivation, nucleation step, lengths, etc.)
        load_case {mustBeMember(load_case, {'P','T','PT'})} = 'P';               % load case 'P': pressure changes, 'T': temperature changes
        load_table table                            % [table] containing time steps, P and T steps (len(y), len(timesteps) for both FW and HW
        stochastic logical = 0;                     % activate stochastic analysis for the single cached member
        diffusion_P logical = 0;                    % activate pressure diffusion
        P_res_mode {mustBeMember(P_res_mode, {'same','different'})} = 'same';                   % base of the reservoir pressure gradient. same = same P at max(depth_HW, depth_FW)
        P0_fault_mode {mustBeMember(P0_fault_mode,{'max','min','mean','FW','HW'})} = 'max';     % [-] assumed initial pressure in fault based on FW and HW pressure. max=max(p_HW, p_FW), etc. 
        P_fault_mode {mustBeMember(P_fault_mode,{'max','min','mean','FW','HW'})} = 'min';       % [-] assumed pressure in fault based on FW and HW pressure during load steps. max=max(p_HW, p_FW), etc. 
        diffusion_T logical = 0;                    % activate pressure diffusion
        aseismic_slip logical = 1                   % compute aseismic slip during nucleation phase
        nucleation_criterion {mustBeMember(nucleation_criterion,{'fixed','UR2D','Day3D','Ruan3D'})} = 'UR2D';   
        nucleation_length_fixed double = 10;           
        save_stress cell = {'all'};                 % indicate which stress to save. 'all', 'none', 'first','last',[step_numbers]
        realizationStale logical = true             % true when faultRealization must be regenerated. Note that it will always be regenerated when executing run.  
        resultsStale logical = true                 % true when stored results no longer match the current configuration
        keepModelObjects logical = false            % true: retain full Pressure/Temperature/FaultStress/FaultSlip objects in *_store
                                                    % after run (for object-method access). Note: raw arrays are always in faultResults
                                                    % regardless; this flag trades extra memory for access to object methods.

    end 

    properties (Constant)
        dx double  = 0;                             % [m] distance from from (for now only on fault allowed)
        n_stochastic {mustBeInteger} = 1;           % retained for compatibility; FaultAnalyzer stores one member only
    end

    properties (Dependent) 
        realizationTable table                      % faultRealization parameters as a flat table for inspection
        nTimes (1,1) double                         % number of modeled load steps
        summary table                               % deprecated alias for faultSummary
    end

    properties (Access = private)
        % Internal result caches. Access via faultResults (arrays) or the
        % Dependent properties (objects, when keepModelObjects = true).
        pressure_store cell = {}                    % {Pressure} when keepModelObjects=true; empty otherwise (data lives in faultResults)
        temperature_store cell = {}                 % {Temperature} when keepModelObjects=true; empty otherwise
        stress_store cell = {}                      % {FaultStress} when keepModelObjects=true; empty otherwise
        slip_store cell = {}                        % always populated after run: {FaultSlip} (keepModelObjects=true) or {slip_meta struct} (lightweight)
        resultsConfiguration = struct()             % configuration snapshot used to produce faultResults
    end

    properties (Dependent, Hidden)
        % Backward-compatibility accessors. Accesss results through
        % getters, and objects to the *_stores
        input_parameters                            % Deprecated. Use faultParameterSpecs instead.
        pressure cell                               % Deprecated. Use faultResults.P. Returns {Pressure obj} if keepModelObjects=true, else struct view.
        temperature cell                            % Deprecated. Use faultResults.T. Returns {Temperature obj} if keepModelObjects=true, else struct view.
        stress cell                                 % Deprecated. Use faultResults.sne/tau. Returns {FaultStress obj} if keepModelObjects=true, else struct view.
        slip cell                                   % Deprecated. Use faultResults.slip. Returns {FaultSlip obj} if keepModelObjects=true, else struct with scalar slip metadata.
        ensemble_members                            % Deprecated. Use faultRealization instead.
        ensemble                                    % Deprecated. Use realizationTable instead.
    end

    methods
        
        function self = FaultAnalyzer(~)
            % PantherInput Load default input parameters
            self.faultParameterSpecs = FaultParameterList(); 
            % delay heavy load_table initialization when performing bulk creation
            % if create_ensemble
            self.load_table = initialize_load_table();
            % self.generateRealization();
        end

        function self = run(self)
            % run Compute stress, slip and nucleation for this fault.
            % Delegates to the three-step extract/compute/apply pattern so
            % MultiFaultAnalyzer can run the compute step in a parfor
            % without broadcasting the full FaultAnalyzer object.
            inputs  = self.extractInputs();
            results = FaultAnalyzer.computeStressAndNucleation(inputs);
            self    = self.applyResults(results);
        end

        function inputs = extractInputs(self)
            % extractInputs Prepare all data needed for run and extract
            % as a plain struct.  Pre-computes Pressure and Temperature so
            % the heavy compute step (computeStressAndNucleation) has no
            % dependency on the FaultAnalyzer object.
            self = self.generateRealization();

            dip       = self.getInputParameter('dip');
            f_s       = self.getDepthDependentInputParameter('f_s');
            f_d       = self.getDepthDependentInputParameter('f_d');
            d_c       = self.getDepthDependentInputParameter('d_c');
            cohesion  = self.getDepthDependentInputParameter('cohesion');

            pressure_obj    = Pressure(self);
            temperature_obj = Temperature(self, 'min');

            inputs = struct();
            inputs.faultRealization     = self.faultRealization;
            inputs.y                   = self.y;
            inputs.dx                  = self.dx;
            inputs.load_case           = self.load_case;
            inputs.nFaultCells         = self.faultLen;
            inputs.nTimeSteps          = self.nTimes;
            % Pre-computed pressure arrays
            inputs.dP_HW = pressure_obj.get_dP_HW();
            inputs.dP_FW = pressure_obj.get_dP_FW();
            inputs.P     = pressure_obj.P;
            inputs.P0    = pressure_obj.P0;
            inputs.dP    = pressure_obj.dP;
            % Pre-computed temperature arrays
            inputs.dT_HW = temperature_obj.get_dT_HW();
            inputs.dT_FW = temperature_obj.get_dT_FW();
            inputs.T     = temperature_obj.T;
            inputs.T0    = temperature_obj.T0;
            inputs.dT    = temperature_obj.dT;
            % Pre-compute Green's functions (geometry only — no time dependence)
            % so workers receive ready-made GF rather than recomputing it.
            [vary_P, vary_T] = FaultStressChange.variableWithDepthGeometryConstant( ....
                inputs.dP_HW, inputs.dP_FW, inputs.dT_HW, inputs.dT_FW);
            vary_dip = FaultStressChange.variableWithDepth(inputs.faultRealization);
            lc = self.load_case;
            vary_PT = (contains(lc,'P') && vary_P) || (contains(lc,'T') && vary_T);
            inputs.GF = GreensFunctions.initialize( ....
                inputs.faultRealization, inputs.y, inputs.dx, vary_PT, vary_dip);
            % Friction / nucleation parameters
            inputs.f_s                      = f_s;
            inputs.f_d                      = f_d;
            inputs.d_c                      = d_c;
            inputs.cohesion                 = cohesion;
            inputs.dip                      = dip;
            inputs.aseismic_slip            = self.aseismic_slip;
            inputs.nucleation_criterion     = self.nucleation_criterion;
            inputs.nucleation_length_fixed  = self.nucleation_length_fixed;
            inputs.keepModelObjects         = self.keepModelObjects;
            % Store model objects only when explicitly requested
            if self.keepModelObjects
                inputs.pressure_obj    = pressure_obj;
                inputs.temperature_obj = temperature_obj;
            end
        end

        function self = applyResults(self, results)
            % applyResults Store computeStressAndNucleation output back
            % into this FaultAnalyzer and refresh the fault summary.
            % Ensure faultRealization is populated (may be empty if this
            % fault object was never run directly via run()).
            if isempty(self.faultRealization) || self.realizationStale
                self = self.generateRealization();
            end
            self.faultResults = results.faultResults;
            self.slip_store   = {results.slip_meta};
            if results.keepModelObjects
                self.pressure_store    = {results.pressure_obj};
                self.temperature_store = {results.temperature_obj};
                self.stress_store      = {results.stress_obj};
                self.slip_store        = {results.slip_obj};
            else
                self.pressure_store    = {};
                self.temperature_store = {};
                self.stress_store      = {};
            end
            self = self.makeResultSummary();
            self.resultsConfiguration = self.captureResultsConfiguration();
            self.resultsStale = false;
        end

        function self = markRealizationStale(self)
            % markRealizationStale [logical] Mark the current faultRealization stale
            % when input parameters have been updated after generating it.
            % faultRealization can be regenerated using generateRealization
            % faultRealization is always regenerated before running
            self.realizationStale = true;
            self.resultsStale = true;
        end

        function self = markResultsStale(self)
            % markResultsStale [logical]] Mark stored results as outdated without deleting them.
            self.resultsStale = true;
        end
       

        function self = setInputParameter(self, parameterName, parameterValues, parameterType)
            if nargin < 4
                parameterType = 'value';
            end
            parameterName = self.validateInputParameterName(parameterName);
            % defensive: when assigning to 'value' ensure a scalar is provided
            if strcmp(parameterType, 'value')
                if ~(isnumeric(parameterValues) && isscalar(parameterValues))
                    error('Assigning to input parameter ''%s'' value must be a numeric scalar', parameterName);
                end
            end
            p = self.faultParameterSpecs.(parameterName);
            p.(parameterType) = parameterValues;
            self.faultParameterSpecs.(parameterName) = p;
            self.realizationStale = true;
            self.resultsStale = true;
        end

        function self = setDepthDependentInputParameter(self, parameterName, parameterValues)
            parameterName = self.validateInputParameterName(parameterName);
            if ~isvector(parameterValues)
                error('parameterValues must be a vector');
            end
            p = self.faultParameterSpecs.(parameterName);
            p.uniform_with_depth = 0;
            p.value_with_depth = parameterValues;
            self.faultParameterSpecs.(parameterName) = p;
            self.realizationStale = true;
            self.resultsStale = true;
        end

        function self = deactivateDepthDependentInputParameter(self, parameterName)
            % deactivateDepthDependentInputParameter Switches an input
            % parameter back to uniform-with-depth mode.
            % Input:
            % parameterName: one of the allowable depth dependent input
            % parameters, e.g. dip
            parameterName = self.validateInputParameterName(parameterName);
            p = self.faultParameterSpecs.(parameterName);
            p.uniform_with_depth = 1;
            self.faultParameterSpecs.(parameterName) = p;
            self.realizationStale = true;
            self.resultsStale = true;
        end

        function self = generateRealization(self)
            % generateRealization  Build the single FaultRealization for this fault.
            % Call this (or run()) before accessing faultRealization.
            self.faultRealization = FaultRealization(self.faultParameterSpecs, self.stochastic);
            self.realizationStale = false;
            self.resultsStale = true;
        end

        function self = generate_ensemble(self) %#ok<MANU>
            % generate_ensemble  Deprecated. Use generateRealization() instead.
            warning('FaultAnalyzer:deprecated', ...
                'generate_ensemble() is deprecated. Use generateRealization() instead.');
            self = self.generateRealization();
        end

        function realizationTable = toRealizationTable(self)
            % toRealizationTable  Return faultRealization parameters as a flat table.
            if isempty(self.faultRealization) || self.realizationStale
                self = self.generateRealization();
            end
            realizationTable = self.faultRealization.to_table();
        end

        function self = makeResultSummary(self)
            %warning('off');
            % reactivation: [boolean] 1 if reactivation detected during any time step, 0 if not
            % reactivation_load_step: [index] index in time array at which
            % reactivation was detected
            % reactivation_dP: [MPa] corresponding pressure change at which
            % reactivation occurred
            % reactivation_dT: [deg] corresponding temperature change at
            % which reactivation occurred
            % nucleation: [boolean] 1 if nucleation detected during any time step, 0 if not
            % nucleation_load_step: [index] index in time array at which
            % nucleation was detected
            % nucleation_dP: [MPa] corresponding pressure change at which
            % nucleation occurred
            % nucleation_dT: [deg] corresponding temperature change at
            % which nucleation occurred
            column_names = {'reactivation', 'reactivation_load_step', 'nucleation', ...
                'nucleation_load_step', 'nucleation_length', 'nucleation_zone_ymid', ...
                'max_slip_length'};
            num_rows = 1;
            self.faultSummary = table(nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                'VariableNames', column_names);
            for i = 1 : length(self.stress)
                self.faultSummary.reactivation(i) = self.slip{i}.reactivation;
                self.faultSummary.reactivation_load_step(i) = self.slip{i}.reactivation_load_step;
                self.faultSummary.nucleation(i) = self.slip{i}.nucleation;
                self.faultSummary.nucleation_load_step(i) = self.slip{i}.nucleation_load_step;
                self.faultSummary.nucleation_length(i) = self.slip{i}.nucleation_length;
                self.faultSummary.nucleation_zone_ymid(i) = self.slip{i}.nucleation_zone_ymid;
                self.faultSummary.max_slip_length(i) = self.slip{i}.max_slip_length;
            end
            warning('on'); 
        end

        function self = make_result_summary(self)
            % make_result_summary Deprecated alias for makeResultSummary.
            warning('FaultAnalyzer:deprecated', ...
                'make_result_summary is deprecated. Use makeResultSummary instead.');
            self = self.makeResultSummary();
        end
        
        function [geom_table] = getRealizationGeometries(self)
            % getRealizationGeometries  Returns geometric indicators for the
            % fault realization as a table.
            % Output:
            %   geom_table - table with geometry, reservoir indices, etc.
            if ~isempty(self.realizationTable)
                input_geometries = {'depth_mid','thick','throw','width_FW', 'width_HW', 'dip'};
                input_table = self.realizationTable;
                geom_table = input_table(:, input_geometries);
                y = self.y;
                for i = 1 : 1
                    geom_table.y_abs{i} = y + geom_table.depth_mid(i);
                    geom_table.L{i} = self.faultRealization.get_along_fault_length(y);
                    geom_table.y_FW_top(i) = self.faultRealization.y_FW_top();
                    geom_table.y_FW_base(i) = self.faultRealization.y_FW_base();
                    geom_table.y_HW_top(i) = self.faultRealization.y_HW_top();
                    geom_table.y_HW_base(i) = self.faultRealization.y_HW_base();
                    geom_table.i_FW_top(i) = self.faultRealization.i_FW_top(y);
                    geom_table.i_FW_base(i) = self.faultRealization.i_FW_base(y);
                    geom_table.i_HW_top(i) = self.faultRealization.i_HW_top(y);
                    geom_table.i_HW_base(i) = self.faultRealization.i_HW_base(y);
                    geom_table.i_FW{i} = self.faultRealization.i_FW(y);
                    geom_table.i_HW{i} = self.faultRealization.i_HW(y);
                    geom_table.i_reservoir{i} = self.faultRealization.i_reservoir(y);
                end
            else
                warning('FaultAnalyzer:noRealization', 'Realization not yet built, call generateRealization() first.');
            end
        end
        
        function [inputParameterValue] = getInputParameter(self, inputParameterName)
            inputParameterName = self.validateInputParameterName(inputParameterName);
            inputParameterValue = self.faultParameterSpecs.(inputParameterName).value;
        end

        function [depthParameterValues] = getDepthDependentInputParameter(self, inputParameterName)
            inputParameterName = self.validateInputParameterName(inputParameterName);
            p = self.faultParameterSpecs.(inputParameterName);
            if p.uniform_with_depth
                depthParameterValues = ones(size(self.y)) * p.value;
            else
                depthParameterValues = p.value_with_depth;
            end
        end

        function [depthParameterValues] = getDepthDependentInput(self, inputParameterName)
            % Convenience alias for getDepthDependentInputParameter
            depthParameterValues = self.getDepthDependentInputParameter(inputParameterName);
        end

        function absoluteDepth = getDepth(self)
            % getDepth Returns absolute depth values using y and depth_mid.
            % Output:
            %   absoluteDepth - Column vector of absolute depth values.
            absoluteDepth = self.y + self.getInputParameter('depth_mid');
        end

        function [LDip, dLDip] = getLAlongDip(self)
            % getLAlongDip Returns along-fault length and spacing from y and dip.
            y = self.y(:);
            dip = self.getDepthDependentInputParameter('dip');

            if isscalar(dip)
                dip = repmat(dip, size(y));
            else
                dip = dip(:);
            end

            if numel(dip) ~= numel(y)
                error('Length of dip (%d) must match length of y (%d)', numel(dip), numel(y));
            end

            LDip = y ./ sind(dip);

            unique_dL = uniquetol(diff(LDip), 0.001);
            if isscalar(unique_dL)
                dLDip = abs(unique_dL);
            else
                dLDip = abs(diff(LDip));
                dLDip = [dLDip; dLDip(end)];
            end
        end

        function output = getResult(self, resultName, varargin)
            % getResult Return a calculated result or derived output.
            allowStale = self.parseAllowStale(varargin{:});
            allowable_result_names = {'P0','P','dP', 'sne', 'tau', 'sne_reac',...
                'tau_reac','sne_nuc','tau_nuc','T0', 'T','dT','slip','scu', ...
                'dcfs','cfs','dcfs_dt','tau_s','tau_d'}';
            if ~ismember(resultName, allowable_result_names)
                resultnames_cellstring = [append(allowable_result_names, repmat({', '},length(allowable_result_names),1))];
                    error(['result name ', resultName, ' not valid, should be one of ', ...
                     resultnames_cellstring{:}]);
            end
            self.requireRunResults(allowStale);
            if isstruct(self.faultResults) && isfield(self.faultResults, resultName)
                output = self.faultResults.(resultName);
            elseif strcmp(resultName, 'scu')
                output = self.getSCU();
            elseif strcmp(resultName, 'tau_s')
                output = self.getStaticFaultStrength();
            elseif strcmp(resultName, 'tau_d')
                output = self.getDynamicFaultStrength();
            elseif strcmp(resultName, 'cfs')
                output = self.getCFF(self.getInputParameter('f_s'), 0);
            elseif strcmp(resultName, 'dcfs')
                cff = self.getCFF(self.getInputParameter('f_s'), 0);
                output = cff - cff(:,1);
            elseif strcmp(resultName, 'dcfs_dt')
                cff = self.getCFF(self.getInputParameter('f_s'), 0);
                time = self.load_table.time_steps;
                % compute the time derivative (MPa/yr)
                output = gradient(cff, time, 2); 
            else
                error('Output %s is not available in faultResults', resultName);
            end
        end

        function outputAtLoadStep = getResultAtLoadStep(self, resultName, loadStep, varargin)
            % getResultAtLoadStep Return a result at an arbitrary load step.
            %
            % Inputs
            %   resultName - result name accepted by getResult
            %   loadStep - scalar load-step index (can be fractional)
            if ~(ischar(resultName) || (isstring(resultName) && isscalar(resultName)))
                error('resultName must be a string');
            end
            if ~(isnumeric(loadStep) && isscalar(loadStep) && isfinite(loadStep))
                error('loadStep must be a finite numeric scalar');
            end
            if loadStep < 1 || loadStep > self.nTimes
                error('loadStep must be between 1 and nTimes (%d)', self.nTimes);
            end

            output = self.getResult(char(resultName), varargin{:});
            outputAtLoadStep = self.sampleResultAlongDimension(output, 2, loadStep, 1:self.nTimes, 'load step');
        end

        function outputAtY = getResultAtY(self, resultName, yValue, varargin)
            % getResultAtY Return a result interpolated at a model y value.
            if ~(isnumeric(yValue) && isscalar(yValue) && isfinite(yValue))
                error('yValue must be a finite numeric scalar');
            end
            output = self.getResult(char(resultName), varargin{:});
            outputAtY = self.sampleResultAlongDimension(output, 1, yValue, self.y, 'y value');
        end

        function outputAtDepth = getResultAtDepth(self, resultName, depthValue, varargin)
            % getResultAtDepth Return a result interpolated at absolute depth.
            if ~(isnumeric(depthValue) && isscalar(depthValue) && isfinite(depthValue))
                error('depthValue must be a finite numeric scalar');
            end
            output = self.getResult(char(resultName), varargin{:});
            outputAtDepth = self.sampleResultAlongDimension(output, 1, depthValue, self.getDepth(), 'depth value');
        end

        function outputAtLoadStep = get_output_at_load_step(self, resultName, loadStep, varargin)
            % get_output_at_load_step Deprecated alias for getResultAtLoadStep.
            warning('FaultAnalyzer:deprecated', ...
                'get_output_at_load_step is deprecated. Use getResultAtLoadStep instead.');
            outputAtLoadStep = self.getResultAtLoadStep(resultName, loadStep, varargin{:});
        end

        function scu = getSCU(self, f_s, cohesion)
            if nargin < 2 || isempty(f_s)
                f_s = self.getDepthDependentInputParameter('f_s');
            end
            if nargin < 3 || isempty(cohesion)
                cohesion = self.getDepthDependentInputParameter('cohesion');
            end
            self.requireRunResults();
            sne = self.faultResults.sne;
            tau = self.faultResults.tau;
            scu = tau ./ (sne .* f_s + cohesion);
        end

        function scu = get_scu(self, f_s, cohesion)
            % Backward-compatible alias for getSCU.
            scu = self.getSCU(f_s, cohesion);
        end

        function tau_s = getStaticFaultStrength(self)
            self.requireRunResults();
            sne = self.faultResults.sne;
            f_s = self.getDepthDependentInputParameter('f_s');
            cohesion = self.getDepthDependentInputParameter('cohesion');
            tau_s = sne .* f_s + cohesion;
        end

        function tau_d = getDynamicFaultStrength(self)
            self.requireRunResults();
            sne = self.faultResults.sne;
            f_d = self.getDepthDependentInputParameter('f_d');
            cohesion = self.getDepthDependentInputParameter('cohesion');
            tau_d = sne .* f_d + cohesion;
        end

        function cff = getCFF(self, mu, cohesion)
            if nargin < 2 || isempty(mu)
                mu = self.getDepthDependentInputParameter('f_s');
            end
            if nargin < 3 || isempty(cohesion)
                cohesion = self.getDepthDependentInputParameter('cohesion');
            end
            self.requireRunResults();
            sne = self.faultResults.sne;
            tau = self.faultResults.tau;
            cff = tau - (sne .* mu + cohesion);
        end

        function cff = get_cff(self, mu, cohesion)
            % Backward-compatible alias for getCFF.
            cff = self.getCFF(mu, cohesion);
        end

        function realizationTable = get.realizationTable(self)
            % Returns an empty table if the realization has not been built
            % yet, to avoid triggering generateRealization() implicitly
            % (e.g. Variable Explorer, parfor broadcast, display).
            % Call generateRealization() or run() first to populate.
            if isempty(self.faultRealization) || self.realizationStale
                realizationTable = table();
            else
                realizationTable = self.faultRealization.to_table();
            end
        end

        % --- Backward-compat getter for deprecated ensemble property ---
        function t = get.ensemble(self)
            warning('FaultAnalyzer:deprecated', ...
                'ensemble is deprecated. Use realizationTable instead.');
            t = self.realizationTable;
        end

        function specs = get.input_parameters(self)
            warning('FaultAnalyzer:deprecated', ...
                'input_parameters is deprecated. Use faultParameterSpecs instead.');
            specs = self.faultParameterSpecs;
        end

        function self = set.input_parameters(self, specs)
            warning('FaultAnalyzer:deprecated', ...
                'input_parameters is deprecated. Use faultParameterSpecs instead.');
            self.faultParameterSpecs = specs;
        end

        % --- Backward-compat getter/setter for deprecated ensemble_members ---
        function em = get.ensemble_members(self)
            if isempty(self.faultRealization)
                em = {};
            else
                em = {self.faultRealization};
            end
        end

        function self = set.ensemble_members(self, val)
            if ~isempty(val)
                self.faultRealization = val{1};
            end
        end

        function pressure = get.pressure(self)
            if ~isempty(self.pressure_store)
                pressure = self.pressure_store;
            elseif ~self.keepModelObjects && isstruct(self.faultResults) && ~isempty(fieldnames(self.faultResults))
                pressure = {struct('P0', self.faultResults.P0, 'P', self.faultResults.P, 'dP', self.faultResults.dP)};
            else
                pressure = {};
            end
        end

        function self = set.pressure(self, pressure)
            self.pressure_store = pressure;
        end

        function temperature = get.temperature(self)
            if ~isempty(self.temperature_store)
                temperature = self.temperature_store;
            elseif ~self.keepModelObjects && isstruct(self.faultResults) && ~isempty(fieldnames(self.faultResults))
                temperature = {struct('T0', self.faultResults.T0, 'T', self.faultResults.T, 'dT', self.faultResults.dT)};
            else
                temperature = {};
            end
        end

        function self = set.temperature(self, temperature)
            self.temperature_store = temperature;
        end

        function stress = get.stress(self)
            if ~isempty(self.stress_store)
                stress = self.stress_store;
            elseif ~self.keepModelObjects && isstruct(self.faultResults) && ~isempty(fieldnames(self.faultResults))
                stress = {struct( ...
                    'sne', self.faultResults.sne, ...
                    'tau', self.faultResults.tau, ...
                    'sne_reac', self.faultResults.sne_reac, ...
                    'tau_reac', self.faultResults.tau_reac, ...
                    'sne_nuc', self.faultResults.sne_nuc, ...
                    'tau_nuc', self.faultResults.tau_nuc, ...
                    'tau_nu', self.faultResults.tau_nu)};
            else
                stress = {};
            end
        end

        function self = set.stress(self, stress)
            self.stress_store = stress;
        end

        function slip = get.slip(self)
            if ~isempty(self.slip_store) && (~self.keepModelObjects || ~isstruct(self.slip_store{1}))
                % In lightweight mode, slip_store keeps scalar metadata.
                meta = self.slip_store{1};
            else
                meta = struct('reactivation', nan, 'reactivation_load_step', nan, 'nucleation', nan, ...
                    'nucleation_load_step', nan, 'nucleation_length', nan, 'nucleation_zone_ymid', nan, 'max_slip_length', nan);
            end

            if ~isempty(self.slip_store) && self.keepModelObjects
                slip = self.slip_store;
                return;
            elseif ~self.keepModelObjects && isstruct(self.faultResults) && ~isempty(fieldnames(self.faultResults))
                meta = struct('reactivation', nan, 'reactivation_load_step', nan, 'nucleation', nan, ...
                    'nucleation_load_step', nan, 'nucleation_length', nan, 'nucleation_zone_ymid', nan, 'max_slip_length', nan);
                if ~isempty(self.slip_store)
                    meta = self.slip_store{1};
                elseif ~isempty(self.faultSummary)
                    meta = struct( ...
                        'reactivation', self.faultSummary.reactivation(1), ...
                        'reactivation_load_step', self.faultSummary.reactivation_load_step(1), ...
                        'nucleation', self.faultSummary.nucleation(1), ...
                        'nucleation_load_step', self.faultSummary.nucleation_load_step(1), ...
                        'nucleation_length', self.faultSummary.nucleation_length(1), ...
                        'nucleation_zone_ymid', self.faultSummary.nucleation_zone_ymid(1), ...
                        'max_slip_length', self.faultSummary.max_slip_length(1));
                end
                slip = {struct( ...
                    'slip', self.faultResults.slip, ...
                    'reactivation', meta.reactivation, ...
                    'reactivation_load_step', meta.reactivation_load_step, ...
                    'nucleation', meta.nucleation, ...
                    'nucleation_load_step', meta.nucleation_load_step, ...
                    'nucleation_length', meta.nucleation_length, ...
                    'nucleation_zone_ymid', meta.nucleation_zone_ymid, ...
                    'max_slip_length', meta.max_slip_length)};
            else
                slip = self.slip_store;
            end
        end

        function self = set.slip(self, slip)
            self.slip_store = slip;
        end

        function nTimes = get.nTimes(self)
            % Single source of truth: number of rows in load_table.
            if ~isempty(self.load_table) && any(strcmp('time_steps', self.load_table.Properties.VariableNames))
                nTimes = height(self.load_table);
            else
                nTimes = 0;
            end
        end

        function summary = get.summary(self)
            % get.summary Backward-compatible alias for faultSummary.
            persistent warned_summary_get
            if isempty(warned_summary_get)
                warned_summary_get = true;
                warning('FaultAnalyzer:DeprecatedSummaryAlias', ...
                    ['FaultAnalyzer.summary is deprecated and will be removed in a future release. ', ...
                    'Use FaultAnalyzer.faultSummary instead.']);
            end
            summary = self.faultSummary;
        end

        function self = set.summary(self, summary)
            % set.summary Backward-compatible alias for faultSummary.
            persistent warned_summary_set
            if isempty(warned_summary_set)
                warned_summary_set = true;
                warning('FaultAnalyzer:DeprecatedSummaryAlias', ...
                    ['Assigning FaultAnalyzer.summary is deprecated and will be removed in a future release. ', ...
                    'Assign FaultAnalyzer.faultSummary instead.']);
            end
            self.faultSummary = summary;
        end

    end

    methods (Access = private)
        function sampledResult = sampleResultAlongDimension(~, result, dimension, query, coordinates, coordinateName)
            if ~isnumeric(result)
                error('Result must be numeric to sample at a %s.', coordinateName);
            end
            coordinates = coordinates(:);
            if isvector(result)
                if dimension == 2
                    result = result(:)';
                else
                    result = result(:);
                end
            end
            if size(result, dimension) ~= numel(coordinates)
                error('Result size along %s dimension (%d) does not match coordinate length (%d).', ...
                    coordinateName, size(result, dimension), numel(coordinates));
            end
            coordinateStep = diff(coordinates);
            if all(coordinateStep < 0)
                coordinates = flipud(coordinates);
                result = flip(result, dimension);
            elseif any(coordinateStep <= 0)
                error('%s coordinates must be strictly monotonic.', coordinateName);
            end
            if query < coordinates(1) || query > coordinates(end)
                error('%s must be between %.6g and %.6g.', coordinateName, coordinates(1), coordinates(end));
            end
            if dimension == 1
                sampledResult = interp1(coordinates, result, query);
            else
                sampledResult = interp1(coordinates, result', query)';
            end
        end

        function requireRunResults(self, allowStale)
            if nargin < 2
                allowStale = false;
            end
            if isempty(self.faultResults) || ~isstruct(self.faultResults) || isempty(fieldnames(self.faultResults))
                error('Run results are not available. Execute FaultAnalyzer.run() first.');
            end
            if self.resultsStale || isempty(fieldnames(self.resultsConfiguration)) || ...
                    ~isequaln(self.resultsConfiguration, self.captureResultsConfiguration())
                self.resultsStale = true;
                if allowStale
                    warning('FaultAnalyzer:StaleResultsAllowed', ...
                        'Returning stale results because AllowStale=true. Run FaultAnalyzer.run() for current results.');
                else
                    error('FaultAnalyzer:StaleResults', ...
                        'Results are stale because model settings changed. Run FaultAnalyzer.run() again.');
                end
            end
        end

        function allowStale = parseAllowStale(~, varargin)
            allowStale = false;
            if mod(numel(varargin), 2) ~= 0
                error('Options must be specified as name-value pairs.');
            end
            for i = 1:2:numel(varargin)
                if ~(ischar(varargin{i}) || (isstring(varargin{i}) && isscalar(varargin{i}))) || ...
                        ~strcmpi(char(varargin{i}), 'AllowStale')
                    error('Unknown option ''%s''. Supported option: AllowStale.', string(varargin{i}));
                end
                value = varargin{i + 1};
                if ~(islogical(value) && isscalar(value))
                    error('AllowStale must be a logical scalar.');
                end
                allowStale = value;
            end
        end

        function configuration = captureResultsConfiguration(self)
            configuration = struct();
            configuration.dy = self.dy;
            configuration.y_extent = self.y_extent;
            configuration.load_case = self.load_case;
            configuration.load_table = self.load_table;
            configuration.stochastic = self.stochastic;
            configuration.diffusion_P = self.diffusion_P;
            configuration.P_res_mode = self.P_res_mode;
            configuration.P0_fault_mode = self.P0_fault_mode;
            configuration.P_fault_mode = self.P_fault_mode;
            configuration.diffusion_T = self.diffusion_T;
            configuration.aseismic_slip = self.aseismic_slip;
            configuration.nucleation_criterion = self.nucleation_criterion;
            configuration.nucleation_length_fixed = self.nucleation_length_fixed;
            configuration.save_stress = self.save_stress;
            configuration.keepModelObjects = self.keepModelObjects;
            parameterNames = properties(self.faultParameterSpecs);
            for i = 1:numel(parameterNames)
                configuration.faultParameterSpecs.(parameterNames{i}) = ...
                    self.faultParameterSpecs.(parameterNames{i});
            end
        end

        function parameterName = validateInputParameterName(self, parameterName)
            % validateInputParameterName Ensures parameter name is text and
            % exists on faultParameterSpecs.
            if ~(ischar(parameterName) || (isstring(parameterName) && isscalar(parameterName)))
                error('parameterName must be a string');
            end
            parameterName = char(parameterName);

            valid_input_parameter_names = properties(self.faultParameterSpecs);
            if ~ismember(parameterName, valid_input_parameter_names)
                validNames = [append(valid_input_parameter_names, repmat({', '}, length(valid_input_parameter_names), 1))];
                error(['input parameter name ', parameterName, ' not valid, should be one of ', validNames{:}]);
            end
        end

    end

    methods (Static)

        function results = computeStressAndNucleation(inputs)
            % computeStressAndNucleation Compute fault stress, slip and
            % nucleation from a plain inputs struct produced by extractInputs.
            %
            % This is a static method with no FaultAnalyzer dependency so
            % it can be called inside a parfor without broadcasting the full
            % object — only the compact inputs struct is sent to each worker.
            y           = inputs.y;
            L           = y ./ sin(inputs.dip * pi / 180);
            nFaultCells = inputs.nFaultCells;
            nTimeSteps  = inputs.nTimeSteps;

            % Initial stress
            initial_stress = InitialStress(y, inputs.faultRealization);

            % Stress changes (uses pre-computed dP/dT arrays and GF from extract_inputs)
            stress_change = FaultStressChange(nFaultCells, nTimeSteps);
            stress_change = stress_change.calc_stress_changes( ...
                inputs.faultRealization, y, inputs.dx, ...
                inputs.dP_HW, inputs.dP_FW, ...
                inputs.dT_HW, inputs.dT_FW, ...
                inputs.load_case, inputs.GF);

            % Total fault stress
            stress_obj = FaultStress(nFaultCells, nTimeSteps);
            stress_obj = stress_obj.compute_fault_stress(initial_stress, stress_change, inputs.P);

            % Aseismic slip and nucleation
            slip_obj = FaultSlip(size(stress_obj.sne, 1), size(stress_obj.sne, 2));
            if inputs.aseismic_slip
                fault_strength = stress_obj.sne .* inputs.f_s + inputs.cohesion;
                [slip_obj, stress_obj.tau] = slip_obj.calculate_fault_slip(L, stress_obj.sne, stress_obj.tau, ...
                    fault_strength, inputs.faultRealization.get_mu_II);
            end
            slip_obj = slip_obj.detect_nucleation(y, L, stress_obj.sne, stress_obj.tau, ...
                inputs.f_s, inputs.f_d, inputs.d_c, inputs.cohesion, ...
                inputs.faultRealization.get_mu_II, ...
                inputs.nucleation_criterion, inputs.nucleation_length_fixed);

            % Reactivation and nucleation stresses
            stress_obj = stress_obj.get_reactivation_stress(slip_obj.reactivation_load_step);
            stress_obj = stress_obj.get_nucleation_stress(slip_obj.nucleation_load_step);

            % Pack results
            results = struct();
            results.keepModelObjects = inputs.keepModelObjects;
            results.faultResults = struct( ...
                'P0', inputs.P0, 'P', inputs.P, 'dP', inputs.dP, ...
                'T0', inputs.T0, 'T', inputs.T, 'dT', inputs.dT, ...
                'sne', stress_obj.sne, 'tau', stress_obj.tau, ...
                'sne_reac', stress_obj.sne_reac, 'tau_reac', stress_obj.tau_reac, ...
                'sne_nuc', stress_obj.sne_nuc, 'tau_nuc', stress_obj.tau_nuc, ...
                'tau_nu', stress_obj.tau_nuc, 'slip', slip_obj.slip);
            results.slip_meta = struct( ...
                'reactivation', slip_obj.reactivation, ...
                'reactivation_load_step', slip_obj.reactivation_load_step, ...
                'nucleation', slip_obj.nucleation, ...
                'nucleation_load_step', slip_obj.nucleation_load_step, ...
                'nucleation_length', slip_obj.nucleation_length, ...
                'nucleation_zone_ymid', slip_obj.nucleation_zone_ymid, ...
                'max_slip_length', slip_obj.max_slip_length);
            if inputs.keepModelObjects
                results.pressure_obj    = inputs.pressure_obj;
                results.temperature_obj = inputs.temperature_obj;
                results.stress_obj      = stress_obj;
                results.slip_obj        = slip_obj;
            end
        end

    end
end





