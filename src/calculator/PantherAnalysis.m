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
        pressure cell
        temperature cell
        stress cell
        slip cell
        faultResults struct
        faultSummary table
    end

    properties (Constant)
        dx double  = 0;                             % [m] distance from from (for now only on fault allowed)
        n_stochastic {mustBeInteger} = 1;           % retained for compatibility; PantherAnalysis stores one member only
    end

    properties (Dependent) 
        ensemble table                              % ensemble member input translated to a table for convenient use
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
            f_s = self.getInputParameter('f_s');
            f_d = self.getInputParameter('f_d');
            d_c = self.getInputParameter('d_c');
            cohesion = self.getInputParameter('cohesion');
            y = self.y; 
            nFaultCells = self.faultLen;
            nTimeSteps = self.nTimes;
            L = y./sin(dip*pi/180);

            % initial stress
            initial_stress{1} = InitialStress(y, self.ensemble_members{1});
            
            % initialize pressure and temperature as function of time
            self.pressure{1} = Pressure(self);
            self.temperature{1} = Temperature(self, 'min');
            
            % stress changes
            stress_change{1} = FaultStressChange(nFaultCells, nTimeSteps);        % initialize fault stresses for P
            stress_change{1} = stress_change{1}.calc_stress_changes( ...
                self.ensemble_members{1}, y, self.dx, ...
                self.pressure{1}.get_dP_HW(), self.pressure{1}.get_dP_FW(), ...
                self.temperature{1}.get_dT_HW(), self.temperature{1}.get_dT_FW(), ...
                self.load_case);
        
            % stress (initial + change)
            self.stress{1} = FaultStress(nFaultCells, nTimeSteps);
            self.stress{1} = self.stress{1}.compute_fault_stress(initial_stress{1}, stress_change{1}, self.pressure{1}.P);
            
            % fault slip, reactivation, nucleation
            self.slip{1} = FaultSlip(size(self.stress{1}.sne, 1), size(self.stress{1}.sne, 2));
            if self.aseismic_slip
                fault_strength{1} = self.stress{1}.sne.*f_s + cohesion;
                [self.slip{1}, self.stress{1}.tau] = self.slip{1}.calculate_fault_slip(L, self.stress{1}.sne, self.stress{1}.tau, ...
                                                             fault_strength{1}, self.ensemble_members{1}.get_mu_II);
            end
            self.slip{1} = self.slip{1}.detect_nucleation(y, L, self.stress{1}.sne, self.stress{1}.tau, f_s, ...
                                                        f_d, d_c, cohesion,self.ensemble_members{1}.get_mu_II, ...
                                                        self.nucleation_criterion, self.nucleation_length_fixed);
            % clear to save memory
            stress_change{1} = [];
            initial_stress{1}= [];

            % get the fault stressses at onset of reactivation and nucleation 
            self.stress{1} = self.stress{1}.get_reactivation_stress(self.slip{1}.reactivation_load_step);
            self.stress{1} = self.stress{1}.get_nucleation_stress(self.slip{1}.nucleation_load_step);

            % aggregate key outputs in a single struct for convenient access
            self.faultResults = struct( ...
                'P0', self.pressure{1}.P0, ...
                'P', self.pressure{1}.P, ...
                'dP', self.pressure{1}.dP, ...
                'T0', self.temperature{1}.T0, ...
                'T', self.temperature{1}.T, ...
                'dT', self.temperature{1}.dT, ...
                'sne', self.stress{1}.sne, ...
                'tau', self.stress{1}.tau, ...
                'sne_reac', self.stress{1}.sne_reac, ...
                'tau_reac', self.stress{1}.tau_reac, ...
                'sne_nuc', self.stress{1}.sne_nuc, ...
                'tau_nuc', self.stress{1}.tau_nuc, ...
                'tau_nu', self.stress{1}.tau_nuc, ...
                'slip', self.slip{1}.slip);

            if ~self.keepModelObjects
                % Replace heavy result objects with lightweight compatibility
                % structs to lower memory usage while keeping common access
                % patterns (e.g., pressure{1}.P, pressure{1}.dP).
                self.pressure{1} = struct('P0', self.faultResults.P0, 'P', self.faultResults.P, 'dP', self.faultResults.dP);
                self.temperature{1} = struct('T0', self.faultResults.T0, 'T', self.faultResults.T, 'dT', self.faultResults.dT);
                self.stress{1} = struct( ...
                    'sne', self.faultResults.sne, ...
                    'tau', self.faultResults.tau, ...
                    'sne_reac', self.faultResults.sne_reac, ...
                    'tau_reac', self.faultResults.tau_reac, ...
                    'sne_nuc', self.faultResults.sne_nuc, ...
                    'tau_nuc', self.faultResults.tau_nuc, ...
                    'tau_nu', self.faultResults.tau_nu);
                self.slip{1} = struct( ...
                    'slip', self.faultResults.slip, ...
                    'reactivation', self.slip{1}.reactivation, ...
                    'reactivation_load_step', self.slip{1}.reactivation_load_step, ...
                    'nucleation', self.slip{1}.nucleation, ...
                    'nucleation_load_step', self.slip{1}.nucleation_load_step, ...
                    'nucleation_length', self.slip{1}.nucleation_length, ...
                    'nucleation_zone_ymid', self.slip{1}.nucleation_zone_ymid, ...
                    'max_slip_length', self.slip{1}.max_slip_length);
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
            if contains(result_name,'P')
                output = self.pressure{run_nr}.(result_name);
            elseif contains(result_name, 'T')
                output = self.temperature{run_nr}.(result_name);
            elseif strcmp(result_name, 'slip')
                output = self.slip{run_nr}.(result_name);
            elseif strcmp(result_name, 'scu')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.getInputParameter('f_s');
                cohesion = self.get_member_input('cohesion');
                output = tau ./ (sne.*f_s + cohesion);
            elseif strcmp(result_name, 'tau_s')
                sne = self.stress{run_nr}.sne;
                f_s = self.get_member_input('f_s');
                cohesion = self.get_member_input('cohesion');
                output = sne.*f_s + cohesion;
            elseif strcmp(result_name, 'tau_d')
                sne = self.stress{run_nr}.sne;
                f_d = self.get_member_input('f_d');
                cohesion = self.get_member_input('cohesion');
                output = sne.*f_d + cohesion;
            elseif strcmp(result_name, 'cfs')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.get_member_input('f_s');
                output = tau - sne.*f_s;
            elseif strcmp(result_name, 'dcfs')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.getInputParameter('f_s');
                output = (tau - tau(:,1)) - sne(sne - sne(:,1)).*f_s;
            elseif strcmp(result_name, 'dcfs_dt')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.getInputParameter('f_s');
                dcfs = (tau - tau(:,1)) - (sne - sne(:,1)).*f_s;
                cfs = tau - sne.*f_s;
                time = self.load_table.time_steps;
                % compute the time derivative (MPa/yr)
                output = gradient(cfs, time, 2); 
            else
                output = self.stress{run_nr}.(result_name);
            end
        end

        function ensemble = get.ensemble(self)
            ensemble = self.ensemble_to_table();
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
        function [cff_max, cff_ymid] = computeCffRates(~, sne, tau, f_s, cohesion, time_yrs, time_range)
            % computeCffRates Equivalent to FaultStress.get_cff_rates but
            % operates on numeric arrays (works for lightweight structs).
            min_index = time_range(1);
            max_index = time_range(2);
            if length(time_yrs) ~= size(sne, 2)
                error('Time steps do not match number of stress output times');
            end
            cff = tau - (sne .* f_s + cohesion);
            cff = cff(:, min_index:max_index);
            time_yrs = time_yrs(min_index:max_index);
            cff_rate = diff(cff, [], 2) ./ diff(time_yrs)';
            cff_max = max(max(cff_rate));
            i_ymid = ceil(size(sne, 1)/2);
            cff_ymid = mean(cff_rate(i_ymid,:));
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