classdef PantherMember < ModelGeometry
    % intializes an ensemble member - i.e. a single model realization
    % for properties that can be depth-dependent, the data type can be a
    % singel number or an array of length(y)

    properties
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
        hyd_diffusivity     % [m2/s] hydraulic diffusivity
        T_grad              % [K/km] temperature gradient
        T_offset            % [k] offset temperature gradient at y=0
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