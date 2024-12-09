classdef PantherMember 
    % intializes an ensemble member - i.e. a single model realization
    % for properties that can be depth-dependent, the data type can be a
    % singel number or an array of length(y)

    properties
        depth_mid           % [m], negative is down
        dip                 % [deg] degrees from horizontal
        dip_azi             % [deg] degrees from north
        thick               % [m] reservoir thickness
        throw               % [m] vertical fault offset
        width_FW            % [m] width footwall compartment
        width_HW            % [m]width footwall compartment
        young               % [Pa] Young's modulus
        poisson             % [-] Poisson's ratio    
        biot                % [-] Biot coefficient
        therm_exp           % [1/K] thermal expansion coefficient
        sH_dir              % [deg] direction SHmax from north
        sHsh                % [-] ratio SH/Sh
        shsv                % [-] ratio Sh/Sv
        sv_grad             % [MPa/km] Vertical stress gradient
        sv_offset           % [MPa] Offset vertical stress gradient at y=0
        p_grad              % [MPa/km] pressure gradient
        p_offset            % [MPa] offset of pressure gradient at y=0    
        p_over              % [MPa] overpressure in the reservoir
        p_grad_res          % [MPa/km] pressure gradient in reservoir  
        p_factor_HW         % [-] depletion factor hanging wall w.r.t unit or P_step
        p_factor_FW         % [-] depletion factor footwall w.r.t. unit or P_step
        p_factor_fault      % [-] depletion factor footwall w.r.t. unit or P_step
        hyd_diffusivity     % [m2/s] hydraulic diffusivity
        T_grad              % [K/km] temperature gradient
        T_offset            % [k] offset temperature gradient at y=0
        T_factor_HW         % [-] cooling factor hanging wall
        T_factor_FW         % [-] cooling factor footwall
        dT_dy_multiplier    % [deg/m] multiplier dT in reservoir wr.t. reservoir mid. -ve is increasing with depth
        therm_diffusivity 
        f_s
        f_d
        d_c
        cohesion
    end    

    methods
        function self = PantherMember(input_parameters, stochastic_analysis)
            % PantherMember creates ensemble member with default or
            % stochastic input parameters, or depth-dependent input
            % INPUT
            % input_parameters: object with input settings for each input 
            % parameters defined in PantherParameterList
            % stochastic analysis: 0 or 1, specified in PantherInput
            param_names = properties(input_parameters);
            for i = 1 : length(param_names)
                % initialize with default value 
                self.(param_names{i}) = input_parameters.(param_names{i}).value;  
                if ~input_parameters.(param_names{i}).uniform_with_depth
                    self.(param_names{i}) = input_parameters.(param_names{i}).value_with_depth;  
                else
                    % and replace if parameter is a stochastic parameter
                    % if a property is depth-dependent it cannot (yet) be
                    % stochastic
                    if input_parameters.(param_names{i}).stochastic && stochastic_analysis
                        dist = input_parameters.(param_names{i}).distribution;
                        a = input_parameters.(param_names{i}).a;
                        b = input_parameters.(param_names{i}).b;
                        random_sample = random(dist, a, b);
                        self.(param_names{i}) = random_sample;
                    end
                end
            end
        end

        function sampled_value = sample_value(distribution, a, b)
            if strcmp(distribution, 'uniform')
                if b < a
                    error('For uniform distribution upper bound b should be larger than lower bound a');
                end
                sampled_value = random('Uniform', a, b);
            else
                error('Other distributions than uniform not yet implemented in the code');
            end

        end

        function [top_HW_y] = top_HW_y(self)
            % top_HW_y returns depth top hanging wall, relative to mid depth
            % top_HW_y = (self.thick- self.throw)/2;
            top_HW_y = get_top_y(self.thick, self.throw, 'HW');
        end

        function [top_FW_y] = top_FW_y(self)
            % top_FW_y returns depth top footwall, relative to mid depth
            % top_FW_y = (self.thick + self.throw)/2;
            top_FW_y = get_top_y(self.thick, self.throw, 'FW');
        end

        function [base_FW_y] = base_FW_y(self)
            % base_FW_y returns depth base footwall, relative to mid depth
            % base_FW_y = -(self.thick - self.throw)/2;
            base_FW_y = get_base_y(self.thick, self.throw, 'FW');
        end

        function [base_HW_y] = base_HW_y(self)
            % base_HW_y returns depth base hangingwall, relative to mid depth
             base_HW_y = get_base_y(self.thick, self.throw, 'HW');
        end

        function [top_HW_i] = top_HW_i(self, y)
            % top_HW_y returns index of first element in top hanging wall 
            top_HW_i = find(y <= self.top_HW_y, 1, 'first');
        end

        function [base_HW_i] = base_HW_i(self, y)
            % top_FW_i index base footwall 
            base_HW_i = find(y >= self.base_HW_y, 1, 'last');
        end

        function [top_FW_i] = top_FW_i(self, y)
            % top_HW_y returns index of first element in top hanging wall 
            top_FW_i = find(y <= self.top_FW_y, 1, 'first');
        end

        function [base_FW_i] = base_FW_i(self, y)
            % top_FW_i index base footwall 
            base_FW_i = find(y >= self.base_FW_y, 1, 'last');
        end        
       
        function [gamma_P] = get_gamma_P(self)
            % get_gamma_P returns the poro-elastic stress path parameter
            gamma_P = (1 - 2*self.poisson).*self.biot./(1 - self.poisson);
        end

        function [gamma_T] = get_gamma_T(self)
            % get_gamma_T returns the thermo-elastic stress path parameter
            % OUTPUT
            % gamma_T   [MPa/deg] thermo-elastic stress path parameter 
            gamma_T = (self.young * self.therm_exp)./(1 - self.poisson);
        end


        function [mu_II] = get_mu_II(self)
            % shear modulus mode II
            mu_II = self.young ./ (2 * (1 + self.poisson) .* (1 - self.poisson));
        end


    end

end