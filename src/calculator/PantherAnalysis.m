classdef (HandleCompatible) PantherAnalysis < FaultMesh 
    % Object that initializes input, sets run and save settings for Panther, and
    % contains the results

    properties
        input_parameters                            % object containing input parameter settings 
        load_case {mustBeMember(load_case, {'P','T','PT'})} = 'P';               % load case 'P': pressure changes, 'T': temperature changes
        load_table table                            % table containing time steps, P and T steps (len(y), len(timesteps) for both FW and HW
        stochastic logical = 0;                     % activate stochastic analysis for the single cached member
        diffusion_P logical = 0;                    % activate pressure diffusion
        P_res_mode {mustBeMember(P_res_mode, {'same','different'})} = 'same';       % base of the reservoir pressure gradient. same = at max(depth_HW, depth_FW)
        P0_fault_mode {mustBeMember(P0_fault_mode,{'max','min','mean','FW','HW'})} = 'max';     % [-] assumed initial pressure in fault based on FW and HW pressure. max=max(p_HW, p_FW), etc. 
        P_fault_mode {mustBeMember(P_fault_mode,{'max','min','mean','FW','HW'})} = 'min';       % [-] assumed pressure in fault based on FW and HW pressure during load steps. max=max(p_HW, p_FW), etc. 
        diffusion_T logical = 0;                    % activate pressure diffusion
        aseismic_slip logical = 1                   % compute aseismic slip during nucleation phase
        nucleation_criterion {mustBeMember(nucleation_criterion,{'fixed','UR2D','Day3D','Ruan3D'})} = 'UR2D';   
        nucleation_length_fixed double = 10;  
        ensemble_members cell                       % single cached member object stored in a 1x1 cell array
        ensemble_dirty logical = true               % indicate whether the cached member must be regenerated
        parallel logical = 1                        % parallel computing for large number of simulations
        save_stress cell = {'all'};                 % indicate which stress to save. 'all', 'none', 'first','last',[step_numbers]
        suppress_status_output logical = false      % indicate ensemble member calculation 
        keepModelObjects logical = false            % keep full result objects (Pressure/Temperature/Stress/Slip) after run
        faultResults struct
        faultSummary table
    end

    properties (Access = private)
        pressure_store cell = {}
        temperature_store cell = {}
        stress_store cell = {}
        slip_store cell = {}
    end

    properties (Constant)
        dx double  = 0;                             % [m] distance from from (for now only on fault allowed)
        n_stochastic {mustBeInteger} = 1;           % retained for compatibility; PantherAnalysis stores one member only
    end

    properties (Dependent) 
        ensemble table                              % ensemble member input translated to a table for convenient use
        pressure cell                               % compatibility view (or stored objects when keepModelObjects=true)
        temperature cell                            % compatibility view (or stored objects when keepModelObjects=true)
        stress cell                                 % compatibility view (or stored objects when keepModelObjects=true)
        slip cell                                   % compatibility view (or stored objects when keepModelObjects=true)
        nTimes (1,1) double                         % number of modeled load steps
        summary table                               % deprecated alias for faultSummary
    end

    methods
        
        function self = PantherAnalysis(~)
            % PantherInput Load default input parameters
            self.input_parameters = PantherParameterList(); 
            % delay heavy load_table initialization when performing bulk creation
            % if create_ensemble
            self.load_table = initialize_load_table();
            % self.generate_ensemble();
        end

        function self = run(self)
            % run the single PantherAnalysis member
            % refresh the cached member
            self.generate_ensemble();
            % run for one ensemble member (multiple members option will be
            % removed in future release)
            
            % unwrap some input parameters for convenience
            dip = self.getInputParameter('dip');
            f_s = self.getDepthDependentInputParameter('f_s');
            f_d = self.getDepthDependentInputParameter('f_d');
            d_c = self.getDepthDependentInputParameter('d_c');
            cohesion = self.getDepthDependentInputParameter('cohesion');
            y = self.y; 
            nFaultCells = self.faultLen;
            nTimeSteps = self.nTimes;
            L = y./sin(dip*pi/180);

            % initial stress
            initial_stress{1} = InitialStress(y, self.ensemble_members{1});
            
            % initialize pressure and temperature as function of time
            pressure_obj = Pressure(self);
            temperature_obj = Temperature(self, 'min');
            
            % stress changes
            stress_change{1} = FaultStressChange(nFaultCells, nTimeSteps);        % initialize fault stresses for P
            stress_change{1} = stress_change{1}.calc_stress_changes( ...
                self.ensemble_members{1}, y, self.dx, ...
                pressure_obj.get_dP_HW(), pressure_obj.get_dP_FW(), ...
                temperature_obj.get_dT_HW(), temperature_obj.get_dT_FW(), ...
                self.load_case);
        
            % stress (initial + change)
            stress_obj = FaultStress(nFaultCells, nTimeSteps);
            stress_obj = stress_obj.compute_fault_stress(initial_stress{1}, stress_change{1}, pressure_obj.P);
            
            % fault slip, reactivation, nucleation
            slip_obj = FaultSlip(size(stress_obj.sne, 1), size(stress_obj.sne, 2));
            if self.aseismic_slip
                fault_strength{1} = stress_obj.sne.*f_s + cohesion;
                [slip_obj, stress_obj.tau] = slip_obj.calculate_fault_slip(L, stress_obj.sne, stress_obj.tau, ...
                                                             fault_strength{1}, self.ensemble_members{1}.get_mu_II);
            end
            slip_obj = slip_obj.detect_nucleation(y, L, stress_obj.sne, stress_obj.tau, f_s, ...
                                                        f_d, d_c, cohesion,self.ensemble_members{1}.get_mu_II, ...
                                                        self.nucleation_criterion, self.nucleation_length_fixed);
            % clear to save memory
            stress_change{1} = [];
            initial_stress{1}= [];

            % get the fault stressses at onset of reactivation and nucleation 
            stress_obj = stress_obj.get_reactivation_stress(slip_obj.reactivation_load_step);
            stress_obj = stress_obj.get_nucleation_stress(slip_obj.nucleation_load_step);

            % keep model objects when requested; otherwise expose array
            % outputs via dependent views from faultResults.
            if self.keepModelObjects
                self.pressure_store = {pressure_obj};
                self.temperature_store = {temperature_obj};
                self.stress_store = {stress_obj};
                self.slip_store = {slip_obj};
            else
                self.pressure_store = {};
                self.temperature_store = {};
                self.stress_store = {};
            end

            % aggregate key outputs in a single struct for convenient access
            self.faultResults = struct( ...
                'P0', pressure_obj.P0, ...
                'P', pressure_obj.P, ...
                'dP', pressure_obj.dP, ...
                'T0', temperature_obj.T0, ...
                'T', temperature_obj.T, ...
                'dT', temperature_obj.dT, ...
                'sne', stress_obj.sne, ...
                'tau', stress_obj.tau, ...
                'sne_reac', stress_obj.sne_reac, ...
                'tau_reac', stress_obj.tau_reac, ...
                'sne_nuc', stress_obj.sne_nuc, ...
                'tau_nuc', stress_obj.tau_nuc, ...
                'tau_nu', stress_obj.tau_nuc, ...
                'slip', slip_obj.slip);

            if ~self.keepModelObjects
                % Keep only lightweight scalar slip metadata in backing
                % storage; pressure/temperature/stress/slip arrays are
                % exposed via dependent views from faultResults.
                self.slip_store = {struct( ...
                    'reactivation', slip_obj.reactivation, ...
                    'reactivation_load_step', slip_obj.reactivation_load_step, ...
                    'nucleation', slip_obj.nucleation, ...
                    'nucleation_load_step', slip_obj.nucleation_load_step, ...
                    'nucleation_length', slip_obj.nucleation_length, ...
                    'nucleation_zone_ymid', slip_obj.nucleation_zone_ymid, ...
                    'max_slip_length', slip_obj.max_slip_length)};
            end

            % Keep run() behavior consistent with panther(): always
            % refresh faultSummary after computing outputs.
            self = self.make_result_summary();

            % % reduce output
            % self.pressure{1} = self.pressure{1}.reduce_steps(indices_for_saving);
            % stress{1} = stress{1}.reduce_steps(indices_for_saving);
            % self.temperature{1} = self.temperature{1}.reduce_steps(indices_for_saving);
            % self.slip{1} = self.slip{1}.reduce_steps(indices_for_saving);
        end

        function self = mark_ensemble_dirty(self)
            self.ensemble_dirty = true;
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
            p = self.input_parameters.(parameterName);
            p.(parameterType) = parameterValues;
            self.input_parameters.(parameterName) = p;
            self.ensemble_dirty = true;
        end

        function self = setDepthDependentInputParameter(self, parameterName, parameterValues)
            parameterName = self.validateInputParameterName(parameterName);
            if ~isvector(parameterValues)
                error('parameterValues must be a vector');
            end
            p = self.input_parameters.(parameterName);
            p.uniform_with_depth = 0;
            p.value_with_depth = parameterValues;
            self.input_parameters.(parameterName) = p;
            self.ensemble_dirty = true;
        end

        function self = deactivateDepthDependentInputParameter(self, parameterName)
            % deactivateDepthDependentInputParameter Switches an input
            % parameter back to uniform-with-depth mode.
            parameterName = self.validateInputParameterName(parameterName);

            p = self.input_parameters.(parameterName);
            p.uniform_with_depth = 1;
            self.input_parameters.(parameterName) = p;
            self.ensemble_dirty = true;
        end

        function self = generate_ensemble(self)
            % Generate the single cached PantherMember used by this analysis.
            self.ensemble_members = cell(1, 1);
            self.ensemble_members{1,1} = PantherMember(self.input_parameters, self.stochastic);
            self.ensemble_dirty = false;
        end

        function ensemble_table = ensemble_to_table(self)
            % create table of input parameter values for easy inspection
            if isempty(self.ensemble_members) || self.ensemble_dirty
                self.generate_ensemble();
            end
            ensemble_table = table;
            for j = 1 : length(self.ensemble_members)
                if j == 1 
                    ensemble_table = self.ensemble_members{j,1}.to_table();
                else
                    new_row = self.ensemble_members{j,1}.to_table();
                    ensemble_table = [ensemble_table; new_row];
                end
            end
        end

        function self = make_result_summary(self)
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
            num_rows = length(self.ensemble_members);
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
        
        function [geom_table] = get_ensemble_geometries(self)
            % geo_geometries Returns useful geometrical indicators for all
            % ensemble members
            % Input:
            % Output:
            %   geom_table - table (height ensemble)
            if ~isempty(self.ensemble)
                input_geometries = {'depth_mid','thick','throw','width_FW', 'width_HW', 'dip'};
                input_table = self.ensemble;
                geom_table = input_table(:, input_geometries);
                y = self.y;
                for i = 1 : length(self.ensemble_members)
                    geom_table.y_abs{i} = y + geom_table.depth_mid(i);
                    geom_table.L{i} = self.ensemble_members{i}.get_along_fault_length(y);
                    geom_table.y_FW_top(i) = self.ensemble_members{i}.y_FW_top();
                    geom_table.y_FW_base(i) = self.ensemble_members{i}.y_FW_base();
                    geom_table.y_HW_top(i) = self.ensemble_members{i}.y_HW_top();
                    geom_table.y_HW_base(i) = self.ensemble_members{i}.y_HW_base();
                    geom_table.i_FW_top(i) = self.ensemble_members{i}.i_FW_top(y);
                    geom_table.i_FW_base(i) = self.ensemble_members{i}.i_FW_base(y);
                    geom_table.i_HW_top(i) = self.ensemble_members{i}.i_HW_top(y);
                    geom_table.i_HW_base(i) = self.ensemble_members{i}.i_HW_base(y);
                    geom_table.i_FW{i} = self.ensemble_members{i}.i_FW(y);
                    geom_table.i_HW{i} = self.ensemble_members{i}.i_HW(y);
                    geom_table.i_reservoir{i} = self.ensemble_members{i}.i_reservoir(y);
                end
            else
                warning('Ensemble not yet initialized, run generate_ensemble first');
            end
        end
        
        function [inputParameterValue] = getInputParameter(self, inputParameterName)
            inputParameterName = self.validateInputParameterName(inputParameterName);
            inputParameterValue = self.input_parameters.(inputParameterName).value;
        end

        function [depthParameterValues] = getDepthDependentInputParameter(self, inputParameterName)
            inputParameterName = self.validateInputParameterName(inputParameterName);
            p = self.input_parameters.(inputParameterName);
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
        

        function [output] = get_member_output(self, result_name, run_nr)
            % getter function to conveniently retrieve output
            if nargin < 3
                run_nr = 1;
            end
            allowable_result_names = {'P0','P','dP' 'sne', 'tau', 'sne_reac',...
                'tau_reac','sne_nuc','tau_nuc','T0', 'T','dT','slip','scu', ...
                'dcfs','cfs','dcfs_dt','tau_s','tau_d'}';
            if ~ismember(result_name, allowable_result_names)
                resultnames_cellstring = [append(allowable_result_names, repmat({', '},length(allowable_result_names),1))];
                error(['result name ', result_name, ' not valid, should be one of ', ...
                     resultnames_cellstring{:}]);
            end
            if run_nr ~= 1
                error('PantherAnalysis stores a single run. Use run_nr = 1.');
            end
            if isstruct(self.faultResults) && isfield(self.faultResults, result_name)
                output = self.faultResults.(result_name);
            elseif strcmp(result_name, 'scu')
                output = self.getSCU(run_nr);
            elseif strcmp(result_name, 'tau_s')
                output = self.getStaticFaultStrength(run_nr);
            elseif strcmp(result_name, 'tau_d')
                output = self.getDynamicFaultStrength(run_nr);
            elseif strcmp(result_name, 'cfs')
                output = self.getCFF(run_nr, self.getInputParameter('f_s'), 0);
            elseif strcmp(result_name, 'dcfs')
                cff = self.getCFF(run_nr, self.getInputParameter('f_s'), 0);
                output = cff - cff(:,1);
            elseif strcmp(result_name, 'dcfs_dt')
                cff = self.getCFF(run_nr, self.getInputParameter('f_s'), 0);
                time = self.load_table.time_steps;
                % compute the time derivative (MPa/yr)
                output = gradient(cff, time, 2); 
            else
                error('Output %s is not available in faultResults', result_name);
            end
        end

        function output_at_load_step = get_output_at_load_step(self, output_name, load_step)
            % get_output_at_load_step Return any faultResults field at an
            % arbitrary load step index between 1 and nTimes.
            %
            % Inputs
            %   output_name - field name in self.faultResults
            %   load_step   - scalar load-step index (can be fractional)
            if ~(ischar(output_name) || (isstring(output_name) && isscalar(output_name)))
                error('output_name must be a string');
            end
            output_name = char(output_name);

            if isempty(self.faultResults) || ~isstruct(self.faultResults)
                error('faultResults is empty. Run PantherAnalysis.run() first.');
            end
            if ~isfield(self.faultResults, output_name)
                valid_fields = fieldnames(self.faultResults);
                valid_fields = [append(valid_fields, repmat({', '}, length(valid_fields), 1))];
                error(['Requested output ''', output_name, ''' not found in faultResults. Valid fields: ', valid_fields{:}]);
            end

            if ~(isnumeric(load_step) && isscalar(load_step) && isfinite(load_step))
                error('load_step must be a finite numeric scalar');
            end
            if load_step < 1 || load_step > self.nTimes
                error('load_step must be between 1 and nTimes (%d)', self.nTimes);
            end

            output = self.faultResults.(output_name);
            if ~isnumeric(output)
                error('faultResults.%s must be numeric', output_name);
            end

            % Time-dependent outputs are sampled along dimension 2.
            if isvector(output)
                if numel(output) == self.nTimes
                    x_ind = 1:self.nTimes;
                    output_at_load_step = interp1(x_ind, output(:), load_step);
                else
                    output_at_load_step = output;
                end
                return;
            end

            n_cols = size(output, 2);
            if n_cols == self.nTimes
                x_ind = 1:self.nTimes;
                output_at_load_step = interp1(x_ind, output', load_step)';
            elseif n_cols == 1
                output_at_load_step = output;
            else
                error(['faultResults.', output_name, ' has size (*,%d) which does not match nTimes (%d)'], n_cols, self.nTimes);
            end
        end

        function scu = getSCU(self, run_nr, f_s, cohesion)
            if nargin < 2 || isempty(run_nr)
                run_nr = 1;
            end
            if run_nr ~= 1
                error('PantherAnalysis stores a single run. Use run_nr = 1.');
            end
            if nargin < 3 || isempty(f_s)
                f_s = self.getDepthDependentInputParameter('f_s');
            end
            if nargin < 4 || isempty(cohesion)
                cohesion = self.getDepthDependentInputParameter('cohesion');
            end
            self.requireRunResults();
            sne = self.faultResults.sne;
            tau = self.faultResults.tau;
            scu = tau ./ (sne .* f_s + cohesion);
        end

        function scu = get_scu(self, run_nr, f_s, cohesion)
            % Backward-compatible alias for getSCU.
            scu = self.getSCU(run_nr, f_s, cohesion);
        end

        function tau_s = getStaticFaultStrength(self, run_nr)
            if nargin < 2 || isempty(run_nr)
                run_nr = 1;
            end
            if run_nr ~= 1
                error('PantherAnalysis stores a single run. Use run_nr = 1.');
            end
            self.requireRunResults();
            sne = self.faultResults.sne;
            f_s = self.getDepthDependentInputParameter('f_s');
            cohesion = self.getDepthDependentInputParameter('cohesion');
            tau_s = sne .* f_s + cohesion;
        end

        function tau_d = getDynamicFaultStrength(self, run_nr)
            if nargin < 2 || isempty(run_nr)
                run_nr = 1;
            end
            if run_nr ~= 1
                error('PantherAnalysis stores a single run. Use run_nr = 1.');
            end
            self.requireRunResults();
            sne = self.faultResults.sne;
            f_d = self.getDepthDependentInputParameter('f_d');
            cohesion = self.getDepthDependentInputParameter('cohesion');
            tau_d = sne .* f_d + cohesion;
        end

        function cff = getCFF(self, run_nr, mu, cohesion)
            if nargin < 2 || isempty(run_nr)
                run_nr = 1;
            end
            if run_nr ~= 1
                error('PantherAnalysis stores a single run. Use run_nr = 1.');
            end
            if nargin < 3 || isempty(mu)
                mu = self.getDepthDependentInputParameter('f_s');
            end
            if nargin < 4 || isempty(cohesion)
                cohesion = self.getDepthDependentInputParameter('cohesion');
            end
            self.requireRunResults();
            sne = self.faultResults.sne;
            tau = self.faultResults.tau;
            cff = tau - (sne .* mu + cohesion);
        end

        function cff = get_cff(self, run_nr, mu, cohesion)
            % Backward-compatible alias for getCFF.
            cff = self.getCFF(run_nr, mu, cohesion);
        end

        function [cff_max, cff_ymid] = get_cff_rates(self, time_range, run_nr, mu, cohesion)
            if nargin < 2 || isempty(time_range)
                time_range = [1, self.nTimes];
            end
            if nargin < 3 || isempty(run_nr)
                run_nr = 1;
            end
            if nargin < 4 || isempty(mu)
                mu = self.getInputParameter('f_s');
            end
            if nargin < 5 || isempty(cohesion)
                cohesion = self.getInputParameter('cohesion');
            end
            cff = self.getCFF(run_nr, mu, cohesion);
            min_index = time_range(1);
            max_index = time_range(2);
            cff = cff(:, min_index:max_index);
            time_yrs = self.load_table.time_steps(min_index:max_index);
            cff_rate = diff(cff, [], 2) ./ diff(time_yrs)';
            cff_max = max(max(cff_rate));
            i_ymid = ceil(size(self.faultResults.sne, 1)/2);
            cff_ymid = mean(cff_rate(i_ymid,:));
        end

        function ensemble = get.ensemble(self)
            ensemble = self.ensemble_to_table();
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
                warning('PantherAnalysis:DeprecatedSummaryAlias', ...
                    ['PantherAnalysis.summary is deprecated and will be removed in a future release. ', ...
                    'Use PantherAnalysis.faultSummary instead.']);
            end
            summary = self.faultSummary;
        end

        function self = set.summary(self, summary)
            % set.summary Backward-compatible alias for faultSummary.
            persistent warned_summary_set
            if isempty(warned_summary_set)
                warned_summary_set = true;
                warning('PantherAnalysis:DeprecatedSummaryAlias', ...
                    ['Assigning PantherAnalysis.summary is deprecated and will be removed in a future release. ', ...
                    'Assign PantherAnalysis.faultSummary instead.']);
            end
            self.faultSummary = summary;
        end

    end

    methods (Access = private)
        function requireRunResults(self)
            if isempty(self.faultResults) || ~isstruct(self.faultResults) || isempty(fieldnames(self.faultResults))
                error('Run results are not available. Execute PantherAnalysis.run() first.');
            end
        end

        function parameterName = validateInputParameterName(self, parameterName)
            % validateInputParameterName Ensures parameter name is text and
            % exists on input_parameters.
            if ~(ischar(parameterName) || (isstring(parameterName) && isscalar(parameterName)))
                error('parameterName must be a string');
            end
            parameterName = char(parameterName);

            valid_input_parameter_names = properties(self.input_parameters);
            if ~ismember(parameterName, valid_input_parameter_names)
                validNames = [append(valid_input_parameter_names, repmat({', '}, length(valid_input_parameter_names), 1))];
                error(['input parameter name ', parameterName, ' not valid, should be one of ', validNames{:}]);
            end
        end

    end
end