classdef PantherParameterList
    % PantherParameterList loads default input parameters
    % PantherParam(value, name, unit, uniform_with_depth, value_with_depth,
    % stochastic, distribution, a, b)
    
    properties
        depth_mid = PantherParam(-2925, 'Mid reservoir depth','y_{mid}','m', 1, nan, 0, 'uniform', -2500, -3200);   % [m], negative is down
        dip = PantherParam(70, 'Fault dip','\theta','degrees', 1, nan, 0, 'uniform', 55, 75);                     % [deg] degrees from horizontal
        dip_azi = PantherParam(90, 'Dip azimuth', '\psi','degrees', 1, nan, 0, 'uniform', 60, 120);              % [deg] degrees from north
        thick = PantherParam(200, 'Thickness','h','m', 1, nan, 0, 'uniform', 100,300);                       % [m] reservoir thickness
        throw = PantherParam(50, 'Throw','t','m', 1, nan, 0, 'uniform', 25, 75);                         % [m] vertical fault offset
        width_FW = PantherParam(inf, 'Width of footwall block','w_{FW}','m', 1, nan, 0, 'uniform', 1000,5000 );   % [m] width footwall compartment
        width_HW = PantherParam(inf, 'Width of hanging wall block','w_{HW}','m', 1, nan,0, 'uniform', 1000,5000);% [m]width footwall compartment
        young = PantherParam(15e3, 'Young''s modulus','E','MPa', 1, nan, 0, 'uniform', 10e3, 20e3);           % [MPa] Young's modulus
        poisson = PantherParam(0.15, 'Poisson''s ratio','\nu','-',1, nan, 0, 'uniform', 0.1, 0.2);             % [-] Poisson's ratio    
        biot = PantherParam(1, 'Biot coefficient','\alpha','-', 1, nan, 0, 'uniform',0.7, 1.0);                   % [-] Biot coefficient
        therm_exp = PantherParam(1e-5, 'Thermal expansion coefficient','\eta_T','1/K', 1, nan, 0, 'uniform', 0.9e-5, 1.2e-5);  %[1/K] thermal expansion coefficient
        sH_dir = PantherParam(0, 'Direction of \sigma_H','\sigma_{Hdir}','degrees', 1, nan, 0, 'uniform', 0, 45);                   % [deg] direction SHmax from north
        sHsh = PantherParam(1, '\sigma_H / \sigma_h','\sigma_H / \sigma_h','-',1, nan, 0, 'uniform', 1.01, 1.1);               % [-] ratio SH/Sh
        shsv = PantherParam(0.75, '\sigma_h / \sigma_v','\sigma_h / \sigma_v','-',1, nan, 0, 'uniform', 0.7, 0.75);            % [-] ratio Sh/Sv
        sv_grad = PantherParam(22, 'Vertical stress gradient','\Delta\sigma_v/\Deltay','MPa/km', 1, nan, 0, 'uniform', 21.5, 22.5);      % [MPa/km] Vertical stress gradient
        sv_offset = PantherParam(0, 'Vertical stress offset','\sigma_{v y=0}','MPa', 1, nan, 0, 'uniform', 0, 1);         % [MPa] Offset vertical stress gradient at y=0
        P_grad = PantherParam(10.5, 'Pressure gradient','\Delta_P/\Deltay','MPa/km', 1, nan, 0, 'uniform', 10, 11);         % [MPa/km] pressure gradient
        P_offset = PantherParam(0, 'Pressure gradient offset','P_{y=0}','MPa', 1, nan, 0, 'uniform', 0, 0);        % [MPa] offset of pressure gradient at y=0    
        P_over = PantherParam(0, 'Hydrostatic overpressure','P_{over}','MPa', 1, nan, 0, 'uniform', 0, 2);            % [MPa] hydrostatic overpressure in/below reservoir 
        P_grad_res = PantherParam(10.5, 'Pressure gradient reservoir fluid','\Delta_{Pres}/\Deltay','MPa/km', 1, nan, 0, 'uniform',10.5, 10.5); % [MPa/km] pressure gradient in reservoir  
        hyd_diffusivity = PantherParam(1e-6, 'Hydraulic diffusivity','Hyd. diff.','m2/s', 1, nan, 0, 'uniform', 0.8e-6,2e-6);   % [m2/s] hydraulic diffusivity
        T_grad = PantherParam(31, 'Temperature gradient','\DeltaT/\Deltay',[char(176),'C/km'], 1, nan, 0, 'uniform', 30, 32);   % [K/km] temperature gradient
        T_offset = PantherParam(10, 'Temperature gradient offset','T_{y=0}',[char(176),'C'], 1, nan, 0, 'uniform', 10, 10); % [k] offset temperature gradient at y=0
        dT_dy_multiplier =  PantherParam(0, 'Depth dependent dT multiplier','Depth dependent dT multiplier','-', 1, nan, 0, 'uniform', 1, 1);        % [-] multiply dT as function of y - y_mid. -ve is increasing dT with depth
        therm_diffusivity = PantherParam(1e-6, 'Thermal diffusivity','T diff.','m2/s', 1, nan, 0, 'uniform', 1e-6, 5e-6);   % [m2/s] thermal diffusivity
        f_s =  PantherParam(0.6, 'Static friction coefficient','f_s','-', 1, nan, 0, 'uniform', 0.5, 0.6);      % [-] static friction coefficient
        f_d =  PantherParam(0.45, 'Dynamic friction coefficient','f_d','-', 1, nan, 0, 'uniform',0.35,0.49);    % [-] dynamic friction coefficient
        d_c =  PantherParam(0.005, 'Critical slip distance','d_c','m', 1, nan, 0, 'uniform',  0.002, 0.01);     % [-] critical slip distance
        cohesion = PantherParam(0, 'Cohesion','C','MPa', 1, nan, 0, 'uniform', 0, 5);                         % [MPa] cohesion
    end

    methods
        
        function self = PantherParameterList()
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

    end

end