classdef (HandleCompatible) PantherInput < FaultMesh
    % Initializes input, and sets run and save settings for Panther

    properties
        input_parameters                            % object containing input parameter settings 
        P_res_mode {mustBeMember(P_res_mode, {'same','different'})} = 'same';       % base of the reservoir pressure gradient. same = at max(depth_HW, depth_FW)
        P0_fault_mode {mustBeMember(P0_fault_mode,{'max','min','mean','FW','HW'})} = 'max';     % [-] assumed initial pressure in fault. max=max(p_HW, p_FW), etc. 
        P_fault_mode {mustBeMember(P_fault_mode,{'max','min','mean','FW','HW'})} = 'min';
        % dp_fault_mode {mustBeMember(dp_fault_mode,{'max','max_abs','min', 'min_abs','mean','FW','HW'})} = 'min';     % [-] assumed pressure in fault. max=max(dp_HW, dp_FW), etc. 
        diffusion_P logical = 0;                    % activate pressure diffusion
        diffusion_T logical = 0;                    % activate pressure diffusion
        stochastic logical = 0;                     % activate stochastic analysis
        n_stochastic {mustBeInteger} = 1;           % number of stochastic runs
        save_stress cell = {'all'};                 % indicate which stress to save. 'all', 'none', 'first','last',[step_numbers]
        load_case {mustBeMember(load_case, {'P','T','PT'})} = 'P';               % load case 'P': pressure changes, 'T': temperature changes
        load_table table                            % table containing time steps, P and T steps (len(y), len(timesteps) for both FW and HW
        ensemble_generated = 0;                     % toggle specifying whether model ensemble has been generated
        ensemble
        parallel logical = 1                        % parallel computing for large number of simulations
        aseismic_slip logical = 1                   % compute aseismic slip during nucleation phase
        nucleation_criterion {mustBeMember(nucleation_criterion,{'fixed','UR2D','Day3D','Ruan3D'})} = 'UR2D';   
        nucleation_length_fixed double = 10;    
    end

    properties (Constant)
        dx double  = 0;                             % [m] distance from from (for now only on fault allowed)
    end

    properties (Dependent) 
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
                for i = 1 : self.n_stochastic
                    self.ensemble{i,1} = PantherMember(self.input_parameters, 1);
                end
            else
                self.ensemble{1,1} = PantherMember(self.input_parameters, 0);
            end
            self.ensemble_generated = 1;
        end

        function reset_ensemble(self)
            self.ensemble = [];
            self.ensemble_generated = 0;
        end

        function ensemble_table = ensemble_to_table(self)
            % create table of input parameter values
            props = properties(self.input_parameters);
            ensemble_table = table;
            if self.ensemble_generated 
                for i = 1 : length(props)
                    for j = 1 : length(self.ensemble)
                        % if non-uniform with depth (length parameters > 1)
                        % store as a NaN in the output table for now
                        if length(self.ensemble{j,1}.(props{i})) > 1
                            ensemble_table.(props{i})(j) = NaN;
                        else
                            ensemble_table.(props{i})(j) = self.ensemble{j,1}.(props{i});
                        end
                    end
                end
            else
                disp('Ensemble not yet generated, execute (your_instance_name).generate_ensemble');
            end
        end

    end
end