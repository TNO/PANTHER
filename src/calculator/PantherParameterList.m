classdef PantherParameterList < handle
    % PantherParameterList loads default input parameters
    % PantherParam(value, name, unit, uniform_with_depth, value_with_depth,
    % stochastic, distribution, a, b)
    
    properties
        depth_mid
        dip
        dip_azi
        thick
        throw
        width_FW
        width_HW
        young
        poisson
        biot
        therm_exp
        sH_dir
        sHsh
        shsv
        sv_grad
        sv_offset
        P_grad
        P_offset
        P_over
        P_grad_res
        hyd_diffusivity
        T_grad
        T_offset
        dT_dy_multiplier
        therm_diffusivity
        f_s
        f_d
        d_c
        cohesion
    end

    methods
        
        function self = PantherParameterList()
            % Construct fresh PantherParam instances per list to avoid shared handles
            self.depth_mid = PantherParam(-2925, 'Mid reservoir depth','y_{mid}','m', 1, nan, 0, 'uniform', -2500, -3200);
            self.dip = PantherParam(70, 'Fault dip','\theta','degrees', 1, nan, 0, 'uniform', 55, 75);
            self.dip_azi = PantherParam(90, 'Dip azimuth', '\psi','degrees', 1, nan, 0, 'uniform', 60, 120);
            self.thick = PantherParam(200, 'Thickness','h','m', 1, nan, 0, 'uniform', 100,300);
            self.throw = PantherParam(50, 'Throw','t','m', 1, nan, 0, 'uniform', 25, 75);
            self.width_FW = PantherParam(inf, 'Width of footwall block','w_{FW}','m', 1, nan, 0, 'uniform', 1000,5000 );
            self.width_HW = PantherParam(inf, 'Width of hanging wall block','w_{HW}','m', 1, nan,0, 'uniform', 1000,5000);
            self.young = PantherParam(15e3, 'Young''s modulus','E','MPa', 1, nan, 0, 'uniform', 10e3, 20e3);
            self.poisson = PantherParam(0.15, 'Poisson''s ratio','\nu','-',1, nan, 0, 'uniform', 0.1, 0.2);
            self.biot = PantherParam(1, 'Biot coefficient','\alpha','-', 1, nan, 0, 'uniform',0.7, 1.0);
            self.therm_exp = PantherParam(1e-5, 'Thermal expansion coefficient','\eta_T','1/K', 1, nan, 0, 'uniform', 0.9e-5, 1.2e-5);
            self.sH_dir = PantherParam(0, 'Direction of \sigma_H','\sigma_{Hdir}','degrees', 1, nan, 0, 'uniform', 0, 45);
            self.sHsh = PantherParam(1, '\sigma_H / \sigma_h','\sigma_H / \sigma_h','-',1, nan, 0, 'uniform', 1.01, 1.1);
            self.shsv = PantherParam(0.75, '\sigma_h / \sigma_v','\sigma_h / \sigma_v','-',1, nan, 0, 'uniform', 0.7, 0.75);
            self.sv_grad = PantherParam(22, 'Vertical stress gradient','\Delta\sigma_v/\Deltay','MPa/km', 1, nan, 0, 'uniform', 21.5, 22.5);
            self.sv_offset = PantherParam(0, 'Vertical stress offset','\sigma_{v y=0}','MPa', 1, nan, 0, 'uniform', 0, 1);
            self.P_grad = PantherParam(10.5, 'Pressure gradient','\Delta_P/\Deltay','MPa/km', 1, nan, 0, 'uniform', 10, 11);
            self.P_offset = PantherParam(0, 'Pressure gradient offset','P_{y=0}','MPa', 1, nan, 0, 'uniform', 0, 0);
            self.P_over = PantherParam(0, 'Hydrostatic overpressure','P_{over}','MPa', 1, nan, 0, 'uniform', 0, 2);
            self.P_grad_res = PantherParam(10.5, 'Pressure gradient reservoir fluid','\Delta_{Pres}/\Deltay','MPa/km', 1, nan, 0, 'uniform',10.5, 10.5);
            self.hyd_diffusivity = PantherParam(1e-6, 'Hydraulic diffusivity','Hyd. diff.','m2/s', 1, nan, 0, 'uniform', 0.8e-6,2e-6);
            self.T_grad = PantherParam(31, 'Temperature gradient','\DeltaT/\Deltay',[char(176),'C/km'], 1, nan, 0, 'uniform', 30, 32);
            self.T_offset = PantherParam(10, 'Temperature gradient offset','T_{y=0}',[char(176),'C'], 1, nan, 0, 'uniform', 10, 10);
            self.dT_dy_multiplier =  PantherParam(0, 'Depth dependent dT multiplier','Depth dependent dT multiplier','-', 1, nan, 0, 'uniform', 1, 1);
            self.therm_diffusivity = PantherParam(1e-6, 'Thermal diffusivity','T diff.','m2/s', 1, nan, 0, 'uniform', 1e-6, 5e-6);
            self.f_s =  PantherParam(0.6, 'Static friction coefficient','f_s','-', 1, nan, 0, 'uniform', 0.5, 0.6);
            self.f_d =  PantherParam(0.45, 'Dynamic friction coefficient','f_d','-', 1, nan, 0, 'uniform',0.35,0.49);
            self.d_c = PantherParam(0.005, 'Critical slip distance','d_c','m', 1, nan, 0, 'uniform',  0.002, 0.01);
            self.cohesion = PantherParam(0, 'Cohesion','C','MPa', 1, nan, 0, 'uniform', 0, 5);
        end

        function [depth_dependent_properties] = get_depth_dependent_properties(self)
            props = properties(self);
            depth_dependent_properties = {};
            for i = 1 : length(props)
                if self.(props{i}).uniform_with_depth == false
                    depth_dependent_properties{end + 1} = props{i};
                end
            end
        end

        function [stochastic_properties, indices] = get_stochastic_properties(self)
            props = properties(self);
            stochastic_properties = {};
            k = 0;
            for i = 1 : length(props)
                if self.(props{i}).stochastic
                    k = k + 1;
                    indices(k) = i;
                    stochastic_properties{k} = props{i};
                end
            end
        end

        function [property_list] = get_parameter_property(self, property_name)
            valid_property_names = properties(self.depth_mid);
            if ~ismember(property_name, valid_property_names)
                valid_property_names_cellstring = [append(valid_property_names , repmat({', '},length(valid_property_names ),1))]; 
                error(['input parameter name ', property_name, ' not valid, should be one of ', ...
                     valid_property_names_cellstring{:}]);
            end
            parameter_list = properties(self);
            property_list = cell(length(parameter_list), 1);
            for i = 1 : length(parameter_list)
                property_list{i} = self.(parameter_list{i}).(property_name);
            end
        end

    end

end