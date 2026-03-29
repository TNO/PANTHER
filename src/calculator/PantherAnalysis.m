classdef (HandleCompatible) PantherAnalysis < FaultMesh 
    % Object that initializes input, sets run and save settings for Panther, and
    % contains the results

    properties
        input_parameters                            % object containing input parameter settings 
        load_case {mustBeMember(load_case, {'P','T','PT'})} = 'P';               % load case 'P': pressure changes, 'T': temperature changes
        load_table table                            % table containing time steps, P and T steps (len(y), len(timesteps) for both FW and HW
        stochastic logical = 0;                     % activate stochastic analysis
        n_stochastic {mustBeInteger} = 1;           % number of stochastic runs
        diffusion_P logical = 0;                    % activate pressure diffusion
        P_res_mode {mustBeMember(P_res_mode, {'same','different'})} = 'same';       % base of the reservoir pressure gradient. same = at max(depth_HW, depth_FW)
        P0_fault_mode {mustBeMember(P0_fault_mode,{'max','min','mean','FW','HW'})} = 'max';     % [-] assumed initial pressure in fault based on FW and HW pressure. max=max(p_HW, p_FW), etc. 
        P_fault_mode {mustBeMember(P_fault_mode,{'max','min','mean','FW','HW'})} = 'min';       % [-] assumed pressure in fault based on FW and HW pressure during load steps. max=max(p_HW, p_FW), etc. 
        diffusion_T logical = 0;                    % activate pressure diffusion
        aseismic_slip logical = 1                   % compute aseismic slip during nucleation phase
        nucleation_criterion {mustBeMember(nucleation_criterion,{'fixed','UR2D','Day3D','Ruan3D'})} = 'UR2D';   
        nucleation_length_fixed double = 10;  
        ensemble_members cell                       % cell array of ensemble member objects (can be generated per request, but will also be regenerated when running PANTHER)
        parallel logical = 1                        % parallel computing for large number of simulations
        save_stress cell = {'all'};                 % indicate which stress to save. 'all', 'none', 'first','last',[step_numbers]
        suppress_status_output logical = false      % indicate ensemble member calculation 
        pressure cell
        temperature cell
        stress cell
        slip cell
        summary table
    end

    properties (Constant)
        dx double  = 0;                             % [m] distance from from (for now only on fault allowed)
    end

    properties (Dependent) 
        ensemble table                              % ensemble member input translated to a table for convenient use
    end

    methods
        
        function self = PantherAnalysis()
            % PantherInput Load default input parameters
            self.input_parameters = PantherParameterList();
            self.load_table = initialize_load_table();
        end

        function self = generate_ensemble(self)
            % Generates ensemble of n_stochastic members
            % Input parameters that are stochastic are randomly sampled
            % for each ensemble member. 
            if self.stochastic
                self.ensemble_members = cell(self.n_stochastic, 1);
                for i = 1 : self.n_stochastic
                    self.ensemble_members{i,1} = PantherMember(self.input_parameters, 1);
                end
            else
                self.ensemble_members = cell(1, 1);
                self.ensemble_members{1,1} = PantherMember(self.input_parameters, 0);
            end
        end

        function ensemble_table = ensemble_to_table(self)
            % create table of input parameter values for easy inspection
            props = properties(self.input_parameters);
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
            column_names = {'reactivation', 'reactivation_load_step','reactivation_dP',...
                'reactivation_dT', 'nucleation', 'nucleation_load_step', 'nucleation_dP',...
                'nucleation_dT','nucleation_length','nucleation_zone_ymid',...
                'slip_length','cff_max', 'cff_ymid','ini_sne','ini_tau'};
            num_rows = length(self.ensemble_members);
            self.summary = table(nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                'VariableNames', column_names);
            for i = 1 : length(self.stress)
                self.summary.reactivation(i) = self.slip{i}.reactivation;
                self.summary.reactivation_load_step(i) = self.slip{i}.reactivation_load_step;
                self.summary.nucleation(i) = self.slip{i}.nucleation;
                self.summary.nucleation_load_step(i) = self.slip{i}.nucleation_load_step;
                n_steps = linspace(1,length(self.load_table.time_steps),length(self.load_table.time_steps));
                % get the reactivation pressure and temperatures
                if ~isnan(self.slip{i}.reactivation_load_step)
                    if strcmp(self.load_case,'P')
                        self.summary.reactivation_dP(i) = interp1(n_steps, self.load_table.P_steps, self.slip{i}.reactivation_load_step);
                        self.summary.reactivation_dT(i) = nan;
                    elseif strcmp(self.load_case,'T') 
                        self.summary.reactivation_dT(i) = interp1(n_steps, self.load_table.T_steps, self.slip{i}.reactivation_load_step);
                        self.summary.reactivation_dP(i) = nan;
                    elseif strcmp(self.load_case,'PT')
                        self.summary.reactivation_dT(i) = interp1(n_steps, self.load_table.T_steps, self.slip{i}.reactivation_load_step);
                        self.summary.reactivation_dP(i) = nan;
                    end
                else
                        self.summary.reactivation_dP(i) = nan;
                        self.summary.reactivation_dT(i) = nan;
                end
                % get the nucleation pressure and temperatures
                if ~isnan(self.slip{i}.nucleation_load_step)
                    if strcmp(self.load_case,'P') 
                        self.summary.nucleation_dP(i) = interp1(n_steps, self.load_table.P_steps, self.slip{i}.nucleation_load_step);
                        self.summary.nucleation_dT(i) = nan;
                    elseif strcmp(self.load_case,'T') 
                        self.summary.nucleation_dT(i) = interp1(n_steps, self.load_table.T_steps, self.slip{i}.nucleation_load_step);
                        self.summary.nucleation_dP(i) = nan;
                    elseif strcmp(self.load_case,'PT') 
                        self.summary.nucleation_dT(i) = interp1(n_steps, self.load_table.T_steps, self.slip{i}.nucleation_load_step);
                        self.summary.nucleation_dP(i) = nan;
                    end
                else
                    self.summary.nucleation_dP(i) = nan;
                    self.summary.nucleation_dT(i) = nan;
                end
                self.summary.nucleation_length(i) = self.slip{i}.nucleation_length;
                self.summary.nucleation_zone_ymid(i) = self.slip{i}.nucleation_zone_ymid;
                if self.slip{i}.nucleation
                    self.summary.slip_length(i) = self.summary.nucleation_length(i);
                else
                    self.summary.slip_length(i) = self.slip{i}.max_slip_length;
                end
                [self.summary.cff_max(i), self.summary.cff_ymid(i)]  = self.stress{i}.get_cff_rates(self.ensemble_members{i}.f_s, self.ensemble_members{i}.cohesion, ...
                self.load_table.time_steps, [1, height(self.load_table)]);
                i_ymid = ceil(size(self.stress{i}.sne,1)/2);
                self.summary.ini_sne(i) = self.stress{i}.sne(i_ymid, 1);
                self.summary.ini_tau(i) = self.stress{i}.tau(i_ymid, 1);
            end
            warning('on'); 
        end
        
        function [geom_table] = get_ensemble_geometries(self)
            % geo_geometries Returns useful geometrical indicators for all
            % ensemble members
            % Input:
            % Output:
            %   geom_table - table (height ensemble)
            self.generate_ensemble();
            input_geometries = {'depth_mid','thick','throw','width_FW', 'width_HW', 'dip'};
            input_table = self.ensemble;
            geom_table = input_table(:, input_geometries);
            y = self.y;
            for i = 1 : length(self.ensemble_members)
                geom_table.y_abs{i} = y + geom_table.depth_mid(i);
                geom_table.L{i} = self.ensemble_members{i}.get_along_fault_length(y);
                geom_table.y_FW_top(i) = self.ensemble_members{i}.y_FW_top();
                geom_table.y_FW_base(i) = self.ensemble_members{i}.y_FW_base();
                geom_table.y_HW_base(i) = self.ensemble_members{i}.y_HW_top();
                geom_table.y_HW_base(i) = self.ensemble_members{i}.y_HW_base();
                geom_table.i_FW_top(i) = self.ensemble_members{i}.i_FW_top(y);
                geom_table.i_FW_base(i) = self.ensemble_members{i}.i_FW_base(y);
                geom_table.i_HW_top(i) = self.ensemble_members{i}.i_HW_top(y);
                geom_table.i_HW_base(i) = self.ensemble_members{i}.i_HW_base(y);
                geom_table.i_FW{i} = self.ensemble_members{i}.i_FW(y);
                geom_table.i_HW{i} = self.ensemble_members{i}.i_HW(y);
                geom_table.i_reservoir{i} = self.ensemble_members{i}.i_reservoir(y);
            end
        end
        
        function [input] = get_member_input(self, input_parameter_name, run_nr)
            if nargin < 3
                run_nr = 1;
            end
            valid_input_parameter_names = properties(self.input_parameters);
            if ~ismember(input_parameter_name, valid_input_parameter_names)
                valid_input_parameter_names_cellstring = [append(valid_input_parameter_names , repmat({', '},length(valid_input_parameter_names ),1))]; 
                error(['input parameter name ', input_parameter_name, ' not valid, should be one of ', ...
                     valid_input_parameter_names_cellstring{:}]);
            end
            input = self.ensemble_members{run_nr}.(input_parameter_name);
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
                f_s = self.get_member_input('f_s');
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
                output = sne.*f_ + cohesion;
            elseif strcmp(result_name, 'cfs')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.get_member_input('f_s');
                output = tau - sne.*f_s;
            elseif strcmp(result_name, 'dcfs')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.get_member_input('f_s');
                output = (tau - tau(:,1)) - sne(sne - sne(:,1)).*f_s;
            elseif strcmp(result_name, 'dcfs_dt')
                sne = self.stress{run_nr}.sne;
                tau = self.stress{run_nr}.tau;
                f_s = self.get_member_input('f_s');
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

    end
end