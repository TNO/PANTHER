classdef PantherParameterList
    % PantherParameterList loads default input parameters
    % PantherParam(value, name, unit, uniform_with_depth, value_with_depth,
    % stochastic, distribution, a, b)
    
    properties
        depth_mid = PantherParam(-2925, 'Mid reservoir depth','m', 1, nan, 0, 'uniform', -2500, -3200);   % [m], negative is down
        dip = PantherParam(70, 'Fault dip','degrees', 1, nan, 0, 'uniform', 55, 75);                     % [deg] degrees from horizontal
        dip_azi = PantherParam(90, 'Dip azimuth','degrees', 1, nan, 0, 'uniform', 60, 120);              % [deg] degrees from north
        thick = PantherParam(200, 'Thickness','m', 1, nan, 0, 'uniform', 100,300);                       % [m] reservoir thickness
        throw = PantherParam(50, 'Thickness','m', 1, nan, 0, 'uniform', 25, 75);                         % [m] vertical fault offset
        width_FW = PantherParam(inf, 'Width of footwall block','m', 1, nan, 0, 'uniform', 1000,5000 );   % [m] width footwall compartment
        width_HW = PantherParam(inf, 'Width of hanging wall block','m', 1, nan,0, 'uniform', 1000,5000);% [m]width footwall compartment
        young = PantherParam(15e3, 'Young''s modulus','Pa', 1, nan, 0, 'uniform', 10e3, 20e3);           % [MPa] Young's modulus
        poisson = PantherParam(0.15, 'Poisson''s ratio','-',1, nan, 0, 'uniform', 0.1, 0.2);             % [-] Poisson's ratio    
        biot = PantherParam(1, 'Biot coefficient','-', 1, nan, 0, 'uniform',0.7, 1.0);                   % [-] Biot coefficient
        therm_exp = PantherParam(1e-5, 'Thermal expansion coefficient','1/K', 1, nan, 0, 'uniform', 0.9e-5, 1.2e-5);  %[1/K] thermal expansion coefficient
        sH_dir = PantherParam(0, 'Direction of \sigma_H','degrees', 1, nan, 0, 'uniform', 0, 45);                   % [deg] direction SHmax from north
        sHsh = PantherParam(1, '\sigma_H / \sigma_h','-',1, nan, 0, 'uniform', 1.01, 1.1);               % [-] ratio SH/Sh
        shsv = PantherParam(0.75, '\sigma_h / \sigma_v','-',1, nan, 0, 'uniform', 0.7, 0.75);            % [-] ratio Sh/Sv
        sv_grad = PantherParam(22, 'Vertical stress gradient','MPa/km', 1, nan, 0, 'uniform', 21.5, 22.5);      % [MPa/km] Vertical stress gradient
        sv_offset = PantherParam(0, 'Vertical stress offset','MPa', 1, nan, 0, 'uniform', 0, 1);         % [MPa] Offset vertical stress gradient at y=0
        p_grad = PantherParam(10.5, 'Pressure gradient','MPa/km', 1, nan, 0, 'uniform', 10, 11);         % [MPa/km] pressure gradient
        p_offset = PantherParam(0, 'Pressure gradient offset','MPa', 1, nan, 0, 'uniform', 0, 0);        % [MPa] offset of pressure gradient at y=0    
        p_over = PantherParam(0, 'Hydrostatic overpressure','MPa', 1, nan, 0, 'uniform', 0, 2);            % [MPa] hydrostatic overpressure in/below reservoir 
        p_grad_res = PantherParam(10.5, 'Pressure gradient reservoir fluid','MPa/km', 1, nan, 0, 'uniform',10.5, 10.5); % [MPa/km] pressure gradient in reservoir  
        %p_factor_HW = PantherParam(1, 'Depletion factor hanging wall','-', 1, nan, 0, 'uniform',1 ,1);   % [-] depletion factor hanging wall w.r.t unit or P_step
        %p_factor_FW = PantherParam(1, 'Depletion factor footwall','-', 1, nan, 0, 'uniform',1, 1);       % [-] depletion factor footwall w.r.t. unit or P_step
        %p_factor_fault = PantherParam(1, 'Depletion factor fault','-', 1, nan, 0, 'uniform', 1, 1);      % [-] depletion factor footwall w.r.t. unit or P_step
        hyd_diffusivity = PantherParam(1e-6, 'Hydraulic diffusivity','m2/s', 1, nan, 0, 'uniform', 0.8e-6,2e-6);   % [m2/s] hydraulic diffusivity
        T_grad = PantherParam(31, 'Temperature gradient',[char(176),'C/km'], 1, nan, 0, 'uniform', 30, 32);   % [K/km] temperature gradient
        T_offset = PantherParam(10, 'Temperature gradient offset',[char(176),'C'], 1, nan, 0, 'uniform', 10, 10); % [k] offset temperature gradient at y=0
        T_factor_HW =  PantherParam(1, 'Cooling factor hanging wall','-', 1, nan, 0, 'uniform', 1, 1);    % [-] cooling factor hanging wall
        T_factor_FW =  PantherParam(1, 'Cooling factor footwall','-', 1, nan, 0, 'uniform', 1, 1);        % [-] cooling factor footwall
        dT_dy_multiplier =  PantherParam(0, 'Depth dependent dT multiplier','-', 1, nan, 0, 'uniform', 1, 1);        % [-] multiply dT as function of y - y_mid. -ve is increasing dT with depth
        therm_diffusivity = PantherParam(1e-6, 'Thermal diffusivity','m2/s', 1, nan, 0, 'uniform', 1e-6, 5e-6);   % [m2/s] thermal diffusivity
        f_s =  PantherParam(0.6, 'Static friction coefficient','-', 1, nan, 0, 'uniform', 0.5, 0.6);      % [-] static friction coefficient
        f_d =  PantherParam(0.45, 'Dynamic friction coefficient','-', 1, nan, 0, 'uniform',0.35,0.49);    % [-] dynamic friction coefficient
        d_c =  PantherParam(0.005, 'Critical slip distance','m', 1, nan, 0, 'uniform',  0.002, 0.01);     % [-] critical slip distance
        cohesion = PantherParam(0, 'Cohesion','MPa', 1, nan, 0, 'uniform', 0, 5);                         % [MPa] cohesion
    end     

end