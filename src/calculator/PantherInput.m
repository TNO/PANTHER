classdef PantherInput < handle
    % Initializes input, and sets run and save settings for Panther
    % Consider renaming to ModelSettings or ModelInput for clarity

    properties
        input_parameters                            % object containing input parameter settings
        dy {mustBePositive} = 2;                    % [m] y-spacing
        y_extent {mustBePositive} = 500;            % [m] extent up and down from depth over which stresses are calculated
        dx double = 0;     % [m] y-spacing
        p_res_mode {mustBeMember(p_res_mode, {'same','different'})} = 'same';       % base of the reservoir pressure gradient. same = at max(depth_HW, depth_FW)
        p_fault {mustBeMember(p_fault,{'max','min','mean','FW','HW'})} = 'min';     % [-] assumed pressure in fault. max=max(p_HW, p_FW), etc. 
        diffusion_P logical = 0;                    % activate pressure diffusion
        diffusion_T logical = 0;                    % activate pressure diffusion
        stochastic logical = 0;                     % activate stochastic analysis
        n_stochastic {mustBeInteger} = 1;           % number of stochastic runs
        save_stress cell = {'all'};                 % indicate which stress to save. 'all', 'none', 'first','last',[step_numbers]
        load_case = 'P';                            % load case 'P': pressure changes, 'T': temperature changes (TODO: combine)
        load_table table                            % table containing time steps, P and T steps (len(y), len(timesteps) for both FW and HW
        ensemble_generated = 0;                     % toggle specifying whether model ensemble has been generated
        ensemble                                    % ensemble of n_stochastic members 
        parallel logical = 1                        % parallel computing for large number of simulations
        aseismic_slip logical = 1                   % compute aseismic slip during nucleation phase
        nucleation_criterion {mustBeMember(nucleation_criterion,{'fixed','UR2D','Day3D','Ruan3D'})} = 'UR2D';   
        nucleation_length_fixed double = 10;    
    end

    properties (Dependent)
         y
    end

    methods
        
        function self = PantherInput()
            % PantherInput Load default input parameters
            self.input_parameters = PantherParameterList;
            self.load_table = initialize_load_table();
            disp('Initialized run instance');
            % automatically initialize 1 ensemble member?
        end

        function self = generate_ensemble(self)
            % Generates ensemble of n_stochastic members
            % Input parameters that are stochastic are randomly sampled
            % for each ensemble member. TODO check in case none of the
            % parameters are stochastic
            if self.stochastic
                for i = 1 : self.n_stochastic
                    self.ensemble{i,1} = PantherMember(self.input_parameters, 1);
                end
                disp(['Ensemble of ', num2str(self.n_stochastic), ' members generated']);
            else
                self.ensemble{1,1} = PantherMember(self.input_parameters, 0);
                disp(['Ensemble of 1 member generated']);
            end
            self.ensemble_generated = 1;
        end

        function reset_ensemble(self)
            self.ensemble = [];
            self.ensemble_generated = 0;
        end

        function self = generate_ensemble_from_table(self, input_table)
            % turn off stochastic (for now). Input values not in the input
            % table will have a fixed value. No-depth dependency. 
            self.stochastic = 0; 
            % check for matching table field names
            panther_input_names = properties(self.input_parameters);
            table_names = fields(input_table);
            matching_columns = ismember(table_names, panther_input_names);
            n_ensemble = height(input_table);
            column_indices = [];
            if any(matching_columns)
                column_indices = find(matching_columns);
               % for k = 1 : 

                self.ensemble_generated = 1;
                for i = 1 : n_ensemble
                    
                    self.ensemble{i,1} = PantherMember(self.input_parameters, 0);
                    for j = 1 : length(column_indices)
                        var = table_names{column_indices(j)};
                        self.ensemble{i,1}.(var) = input_table.(var)(i);
                    end
                end
            else
                disp('Table column names do not match Panther input parameter names');
            end
            
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

        function a = get.y(self)
            % y initializes depth y relative to depth_mid
            a = get_fault_y(self.dy, self.y_extent);
        end

        % TODO: add csv writer and reader. include example csv
        % TODO: allow for different thickness HW and FW
        % TODO: allow for custom depth-dependent input parameters (in that case, disable stochastic?)
        % functionality to check whether all stress path parameters are
        % depth-independent (nu, biot, dip)
        % interpolate to y. 
    end
end