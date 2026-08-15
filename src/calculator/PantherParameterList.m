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
            % Build the input parameter structures
            self.depth_mid = makePantherParam(-2925, 'Mid reservoir depth','y_{mid}','m', 1, nan, 0, 'uniform', -2500, -3200);
            self.dip = makePantherParam(70, 'Fault dip','\theta','degrees', 1, nan, 0, 'uniform', 55, 75);
            self.dip_azi = makePantherParam(90, 'Dip azimuth', '\psi','degrees', 1, nan, 0, 'uniform', 60, 120);
            self.thick = makePantherParam(200, 'Thickness','h','m', 1, nan, 0, 'uniform', 100,300);
            self.throw = makePantherParam(50, 'Throw','t','m', 1, nan, 0, 'uniform', 25, 75);
            self.width_FW = makePantherParam(inf, 'Width of footwall block','w_{FW}','m', 1, nan, 0, 'uniform', 1000,5000 );
            self.width_HW = makePantherParam(inf, 'Width of hanging wall block','w_{HW}','m', 1, nan,0, 'uniform', 1000,5000);
            self.young = makePantherParam(15e3, 'Young''s modulus','E','MPa', 1, nan, 0, 'uniform', 10e3, 20e3);
            self.poisson = makePantherParam(0.15, 'Poisson''s ratio','\nu','-',1, nan, 0, 'uniform', 0.1, 0.2);
            self.biot = makePantherParam(1, 'Biot coefficient','\alpha','-', 1, nan, 0, 'uniform',0.7, 1.0);
            self.therm_exp = makePantherParam(1e-5, 'Thermal expansion coefficient','\eta_T','1/K', 1, nan, 0, 'uniform', 0.9e-5, 1.2e-5);
            self.sH_dir = makePantherParam(0, 'Direction of \sigma_H','\sigma_{Hdir}','degrees', 1, nan, 0, 'uniform', 0, 45);
            self.sHsh = makePantherParam(1, '\sigma_H / \sigma_h','\sigma_H / \sigma_h','-',1, nan, 0, 'uniform', 1.01, 1.1);
            self.shsv = makePantherParam(0.75, '\sigma_h / \sigma_v','\sigma_h / \sigma_v','-',1, nan, 0, 'uniform', 0.7, 0.75);
            self.sv_grad = makePantherParam(22, 'Vertical stress gradient','\Delta\sigma_v/\Deltay','MPa/km', 1, nan, 0, 'uniform', 21.5, 22.5);
            self.sv_offset = makePantherParam(0, 'Vertical stress offset','\sigma_{v y=0}','MPa', 1, nan, 0, 'uniform', 0, 1);
            self.P_grad = makePantherParam(10.5, 'Pressure gradient','\Delta_P/\Deltay','MPa/km', 1, nan, 0, 'uniform', 10, 11);
            self.P_offset = makePantherParam(0, 'Pressure gradient offset','P_{y=0}','MPa', 1, nan, 0, 'uniform', 0, 0);
            self.P_over = makePantherParam(0, 'Hydrostatic overpressure','P_{over}','MPa', 1, nan, 0, 'uniform', 0, 2);
            self.P_grad_res = makePantherParam(10.5, 'Pressure gradient reservoir fluid','\Delta_{Pres}/\Deltay','MPa/km', 1, nan, 0, 'uniform',10.5, 10.5);
            self.hyd_diffusivity = makePantherParam(1e-6, 'Hydraulic diffusivity','Hyd. diff.','m2/s', 1, nan, 0, 'uniform', 0.8e-6,2e-6);
            self.T_grad = makePantherParam(31, 'Temperature gradient','\DeltaT/\Deltay',[char(176),'C/km'], 1, nan, 0, 'uniform', 30, 32);
            self.T_offset = makePantherParam(10, 'Temperature gradient offset','T_{y=0}',[char(176),'C'], 1, nan, 0, 'uniform', 10, 10);
            self.dT_dy_multiplier =  makePantherParam(0, 'Depth dependent dT multiplier','Depth dependent dT multiplier','-', 1, nan, 0, 'uniform', 1, 1);
            self.therm_diffusivity = makePantherParam(1e-6, 'Thermal diffusivity','T diff.','m2/s', 1, nan, 0, 'uniform', 1e-6, 5e-6);
            self.f_s =  makePantherParam(0.6, 'Static friction coefficient','f_s','-', 1, nan, 0, 'uniform', 0.5, 0.6);
            self.f_d =  makePantherParam(0.45, 'Dynamic friction coefficient','f_d','-', 1, nan, 0, 'uniform',0.35,0.49);
            self.d_c = makePantherParam(0.005, 'Critical slip distance','d_c','m', 1, nan, 0, 'uniform',  0.002, 0.01);
            self.cohesion = makePantherParam(0, 'Cohesion','C','MPa', 1, nan, 0, 'uniform', 0, 5);
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
            % determine valid property names of a PantherParam (struct or object)
            if isstruct(self.depth_mid)
                valid_property_names = fieldnames(self.depth_mid);
            else
                valid_property_names = properties(self.depth_mid);
            end
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

    methods
        function varargout = subsasgn(self, S, value)
            % Intercept nested assignments like obj.param_name.field = value
            % when param_name is one of the input_parameters. Perform a
            % get-modify-set so the change persists in this handle object.
            try
                if ~isempty(S) && isequal(S(1).type, '.') && numel(S) >= 2 && isequal(S(2).type, '.')
                    prop = S(1).subs;
                    props = properties(self);
                    if ismember(prop, props)
                        p = self.(prop);
                        if isstruct(p)
                            % perform local modification and reassign
                            field = S(2).subs;
                            if numel(S) > 2
                                % There is further indexing into the subfield, e.g.
                                % obj.param.field(i:j) = value. Apply the tail of
                                % the substruct to the subfield using builtin subsasgn
                                tail = S(3:end);
                                oldsub = p.(field);
                                newsub = builtin('subsasgn', oldsub, tail, value);
                                p.(field) = newsub;
                            else
                                p.(field) = value;
                            end
                            self.(prop) = p;
                            if nargout == 0
                                return;
                            else
                                [varargout{1:nargout}] = deal();
                                return;
                            end
                        end
                    end
                end
            catch
                % fall through to builtin behavior on error
            end
            if nargout == 0
                builtin('subsasgn', self, S, value);
            else
                [varargout{1:nargout}] = builtin('subsasgn', self, S, value);
            end
        end
    end

end