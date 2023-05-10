function [varin] = prepare_input_for_stress_calculation(varin)
    
    % elastic parameters
    varin.mu = varin.E/(2*(1+varin.nu));                  % [Pa] shear modulus
    varin.g = 9.81;                 % [m/s^2] gravitational constant
    biot = varin.biot;              % [-] Biot's coefficient 
    E = varin.E;                    % [Pa] Young's modulus
    nu = varin.nu;                  % [-] Poisson's ratio
    mu = varin.mu;                    % [Pa] shear modulus
    thermexp = varin.thermexp;      % [1/K] thermal expansion coefficient

    dpres = 1;                      % [Pa] unit pressure or temperature to compute stress change

    if ~isfield(varin,'load_case')
        load_case = 'P';
    else
        load_case = varin.load_case;
    end
    % coefficients for displacement 
    if strcmp(load_case,'P')
        varin.D = (1-2*nu)*biot*dpres / (2*pi*(1-nu)*mu);   % coefficient for displacements, - 
    elseif strcmp(load_case,'T')
        varin.D = E*thermexp*dpres / (2*pi*(1-nu)*mu);      % coefficient for displacements, - 
    end
    varin.C = mu*varin.D;                    % coefficient for stresses, N/m

    % calculate fault geometry inut for Jansen 2019 stress change calculation
    varin.theta = varin.dip_angle*pi/180;             % [rad] fault dip, positive = normal fault 
    thickness = varin.thick;
%     if varin.wR == 0 || varin.wL ==0
%         varin.throw = 0;
%     end
    offset = varin.throw;                   % set to 0 for wR or wL == 0?
    varin.a = 0.5*thickness - 0.5*offset;   % [m] geometrical descriptor Jasne 2019
    varin.b = (thickness + offset)/2;       % [m] reservoir mid depth (from bottom of footwall or top hw)
    varin.D_center = varin.cdepth;  

    varin.dx = 0;                           % [m] horizontal distance from fault plane (x=0: (almost) on the fault)
    if or(isempty(varin.width), ~isfield(varin, 'width'))
        varin.width = inf;                  % if no width is specified, set to infinite
    end
    if isfield(varin,'wL')
        if ~isempty(varin.wL)
            varin.c = varin.wL;         % width of the right compartment
        end
    else
        varin.c = varin.width;
    end
    if isfield(varin,'wR')
        if ~isempty(varin.wR)
            varin.d = varin.wR;         % width of the left compartment
        end
    else
        varin.d = varin.width;
    end

end