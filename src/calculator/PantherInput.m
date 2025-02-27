classdef (HandleCompatible) PantherInput < FaultMesh
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
        
        function self = PantherInput()
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
            column_names = {'reactivation', 'reactivation_load_step','reactivation_dP',...
                'reactivation_dT', 'nucleation', 'nucleation_load_step', 'nucleation_dP',...
                'nucleation_dT','cff_max', 'cff_ymid','ini_sne','ini_tau'};
            num_rows = length(self.ensemble_members);
            self.summary = table(nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),nan(num_rows,1),...
                'VariableNames', column_names);
            for i = 1 : length(self.stress)
                self.summary.reactivation(i) = self.slip{i}.reactivation;
                self.summary.reactivation_load_step(i) = self.slip{i}.reactivation_load_step;
                self.summary.nucleation(i) = self.slip{i}.nucleation;
                self.summary.nucleation_load_step(i) = self.slip{i}.nucleation_load_step;
                n_steps = linspace(1,length(self.load_table.time_steps),length(self.load_table.time_steps));
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
                [self.summary.cff_max(i), self.summary.cff_ymid(i)]  = self.stress{i}.get_cff_rates(self.ensemble_members{i}.f_s, self.ensemble_members{i}.cohesion, ...
                self.load_table.time_steps, [1, height(self.load_table)]);
                i_ymid = ceil(size(self.stress{i}.sne,1)/2);
                self.summary.ini_sne(i) = self.stress{i}.sne(i_ymid, 1);
                self.summary.ini_tau(i) = self.stress{i}.tau(i_ymid, 1);
            end
            warning('on'); 
        end


        function ensemble = get.ensemble(self)
            ensemble = self.ensemble_to_table();
        end

    end
end