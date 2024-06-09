classdef MultiFaultCalculator
    % object to perform calculations on the multiple faults (called pillars) 
    % for each pillar different input and different run settings can be
    % specified
    
    properties
        pillars cell        % cell array of PANTHER input objects (1 ensemble per entry, can be modified to multiple) 
        pillar_info table   % table with pillar custom meta_data (e.g. name, coordinates)
        pillar_results cell
        result_summary table
        run_done logical
    end

    properties (Dependent)
        n_pillars double
    end

    methods

        function self = MultiFaultCalculator(n_pillars)
            % construct the class with n_pillars, assign ID in the
            % information table
            self.pillars = cell(n_pillars, 1);
            self.pillar_info = table([1:n_pillars]','VariableNames',{'ID'});
            % initialize the default PANTHER input for each pillar
            for i = 1 : length(self.pillars)
                self.pillars{i} = PantherInput();
            end
        end


        function self = run(self)
            % run simulation for all the fault pillars
            tic
            all_pillars = self.pillars;
            n = self.n_pillars;
            self.pillar_results = cell(n, 1);
            parfor i = 1 : n
                results{i,1} = panther(all_pillars{i});
                disp([num2str(i),'/', num2str(n)]);
            end
            self.pillar_results = results;
            toc
            self.run_done = true;
            self.result_summary = self.get_results_summary();
        end
        

        function self = add_pillar_info_as_table(self, info_table_to_be_added)
            % adds meta data into the pillar_info table
            % INPUT
            % info_table_to_be_added    table of height n_faults
            if height(info_table_to_be_added) ~= height(self.pillar_info)
                disp(['Cant append fault info, table size does not match. Height should be ',num2str(height(self.pillar_info)) ]);
            else
                self.pillar_info = [self.pillar_info, info_table_to_be_added];
            end
        end


        function self = set_depth_dependent_input_parameter(self, parameter_name, values)
            % sets numeric input Panther input parameters
            % INPUT
            % values    cell array length(n_pillars), containing arrays of
            % doubles of length(y)
            if ~iscell(values)
                error('Input depth dependent variable in a cell array of length n_pillars');
            end
            if self.is_valid_input_parameter_name(parameter_name) && (length(values) == length(self.pillars))
                for i = 1 : length(self.pillars)
                    if length(values{i}) == length(self.pillars{i}.y)
                        self.pillars{i}.input_parameters.(parameter_name).uniform_with_depth = 0;
                        self.pillars{i}.input_parameters.(parameter_name).value_with_depth = values{i};
                    else
                        disp(['depth dependent variable could not be set, size not equal to y, length is ', num2str(self.pillars{i}.y)]);
                    end
                end
            else
                disp(['variable ', parameter_name,' not assigned, check that length of input values equals number of pillars']);
            end
        end


        function self = set_input_parameter(self, parameter_name, values, property)
            % sets numeric input Panther input parameters
            % INPUT
            % values    array of doubles
            % property  'value', 'a', or 'b' 
            if nargin < 4
                property = 'value';
            end
            if self.is_valid_input_parameter_name(parameter_name) && (length(values) == length(self.pillars))
                for i = 1 : length(self.pillars)
                    self.pillars{i}.input_parameters.(parameter_name).(property) = values(i);
                end
            else
                disp(['variable ', parameter_name,' not assigned, check that length of input values equals number of pillars']);
            end
        end


        function [input_values] = get_input_value(self, parameter_name, input_property)
            % retrieve input parameters from the fault pillars
            % INPUT
            % parameter_name    string. parameter name, e.g. 'dip'
            if nargin < 4
                input_property = 'value';
            end
            input_values = zeros(length(self.pillars), 1);
            if self.is_valid_input_parameter_name(parameter_name)
                 for i = 1 : length(self.pillars)
                    input_values(i) = self.pillars{i}.input_parameters.(parameter_name).(input_property);
                end
            end
        end
        
        function self = set_run_setting(self, setting_name, setting_value)
            % specify run settings per pillar. (exc load table, y, input)
            % INPUT
            % setting name
            % setting value     cell array, or array of floats or single
            % cell/single float
            [is_valid, value_type] = self.is_valid_setting_name(setting_name);
            if is_valid
                for i = 1 : length(self.pillars)
                    valid_value = 1;
                    if iscell(setting_value) && (length(setting_value) == length(self.pillars))
                        assign_value = setting_value{i};
                    elseif iscell(setting_value) && (length(setting_value) == 1)
                        assign_value = setting_value{1};
                    elseif isfloat(setting_value) && (length(setting_value) == 1)
                        assign_value = setting_value(1);    
                    elseif isfloat(setting_value) && (length(setting_value) == self.n_pillars)
                        assign_value = setting_value(i);  
                    else
                        valid_value = 0;
                    end
                    if valid_value & strcmp(value_type,class(assign_value))
                        self.pillars{i}.(setting_name) = assign_value;
                    end
                end
                if ~valid_value || ~strcmp(value_type,class(assign_value))
                    disp('Specified setting type does not seem the right type or dimension, check');
                end
            end

        end


        function [nuc_load_step] = get_minimum_nucleation_load_step(self)
            if self.run_done
                nuc_load_step = min(self.result_summary.nucleation_load_step);
            else
                nuc_load_step = nan;
            end
        end


        function self = overwrite_nucleation_stress(self, new_nucleation_load_step)
            for i = 1 : length(self.pillar_results)
                nuc = new_nucleation_load_step;
                self.pillar_results{i}.stress{1} = self.pillar_results{i}.stress{1}.get_nucleation_stress(nuc);
            end
        end

        function self = reduce_output(self, time_step_indices)
            % return output only at give time step indices
            
            for i = 1 : length(self.pillar_results)
%                 if max(time_step_indices) < size(self.pillar_results{i}.stress{1}.sne, 2)  & ...
%                 (min(time_step_indices) > 1)
                    self.pillar_results{i}.stress{1} = self.pillar_results{i}.stress{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.temperature{1} = self.pillar_results{i}.temperature{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.pressure{1} = self.pillar_results{i}.pressure{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.slip{1} = self.pillar_results{i}.slip{1}.reduce_steps(time_step_indices);
                    self.pillar_results{i}.load_table = self.pillar_results{i}.load_table(time_step_indices,:);
                end
%             end
        end

        function self = add_info_from_closest_point(self, X, Y, Z_value, X_column, Y_column, new_column)
            % add meta data info based on nearest point 
            % INPUT
            % X             x-coordinates
            % Y             y-coordinates
            % Z_value       spatial value
            % X_column      string, column name of the X_coordinate in pillar_info
            % Y_column      string, column name of the X_coordinate in pillar_info
            % new_column    string, column name of Z_value
            % TODO add some checks
            i5 = floor(self.n_pillars/20);
            reverseStr = '';
            for i = 1 : length(self.pillars)
                xq = self.pillar_info.(X_column)(i);
                yq = self.pillar_info.(Y_column)(i);
                [closest_index] = self.nearest_rectangular_grid_coordinate(X,Y,xq, yq);
                if ~isempty(closest_index)
                    self.pillar_info.(new_column)(i) = Z_value(closest_index(1));
                    %self.display_progress(i, self.n_pillars, 20, ...
                    %    ['Progress assigning spatial data ', new_column,' : ']);
                    if rem(i, i5) == 0
                        percentDone = 100*i / self.n_pillars;
                        msg = sprintf('Progress assigning spatial data : %3.1f', percentDone);
                        fprintf([reverseStr, msg]);
                        reverseStr = repmat(sprintf('\b'), 1, length(msg));
                    end
                end
            end
            sprintf(newline);
        end
        
        function [summary] = get_results_summary(self)
            if self.run_done
                summary = self.pillar_results{1}.summary;
                % concatenate summary tables of individual fault pillars
                for i = 1 : self.n_pillars - 1
                    summary = [summary; self.pillar_results{i}.summary];
                end
            else
                disp('Run not yet exectued, empty summary');
                summary = [];
            end
        end
      
        function [i_min] = nearest_rectangular_grid_coordinate(~, X_grid, Y_grid, X_query, Y_query)
            % only works for rectangular grids
            [~, i_min_x] = min(abs(X_query - X_grid));
            [~, i_min_y] = min(abs(Y_query - Y_grid));
            i_min = find((X_grid == X_grid(i_min_x)) & (Y_grid == Y_grid(i_min_y)));
        end

        function [valid_name] = is_valid_input_parameter_name(self, submitted_name)
            % validate whether specified input parameter name is valid
            valid_field_names = fields(self.pillars{1}.input_parameters);
            if ismember(submitted_name, valid_field_names)
                valid_name = true;
            else
                valid_name = false;
                fields_cellstring = [append(valid_field_names, repmat({', '},length(valid_field_names),1))];
                disp(['Check length, or input parameter name not valid, should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        function [valid_name, value_type] = is_valid_setting_name(self, submitted_name)
            % check if run setting name is valid
            valid_setting_names = fields(self.pillars{1});
            if ismember(submitted_name, valid_setting_names) & ~ismember(submitted_name,{'input_parameters','load_table','y','ensemble'})
                valid_name = true;
                if ismember(submitted_name,{'p_res_mode','p_fault','load_case','nucleation_criterion'})
                    value_type = 'char';
                else
                    value_type = 'double';
                end
            else
                valid_name = false;
                value_type = 'double';
                fields_cellstring = [append(valid_setting_names, repmat({', '},length(valid_setting_names),1))];
                disp(['Run setting name not valid, should be one of the following: ',...
                     [fields_cellstring{:}]]);
            end
        end

        
%         function display_progress(~, iteration, total_number, steps, msg_prefix)
%             if nargin < 5
%                 msg_prefix = 'Percent done ';
%             end
%             i_step = floor(total_number/steps);
%             reverseStr = '';
%             if rem(iteration, i_step) == 0
%                 percentDone = 100*(iteration / total_number);
%                 % msg = sprintf([msg_prefix,  num2str(percentDone)]);
%                 msg = sprintf('Progress reading porosity: %3.1f', percentDone);
%                 fprintf([reverseStr, msg]);
%                 reverseStr = repmat(sprintf('\b'), 1, length(msg));
%             end
%         end
        
        function num = get.n_pillars(self)
            % number of faults in the object
            num = length(self.pillars);
        end

    end

end
