classdef InitialStress
    % InitialStress computes the initial, pre-production total normal stress 
    % and shear stress, based on the depth and specified stress gradients, 
    % ratios and fault geometry

    properties
        sn0 double       % [MPa] initial total normal stress
        tau0 double      % [MPa] initial shear stress
    end

    methods
        function self = InitialStress(y, params)
            % Initializes initial stress
            % INPUT 
            % y:        depth with respect to 0,0
            % params:   input parameters for the model member
            y_abs = y + params.depth_mid;
            Sv = -(y_abs/1000) .* params.sv_grad + params.sv_offset;    % [MPa] vertical stress
            Sh = params.shsv .* Sv;         % [MPa] total minimum horizontal stress
            SH = Sh.*params.sHsh;           % [MPa] total maximum horizontal stress
            SHdir = params.sH_dir;          % [deg] direction of maximum horizontal stress
            dip = params.dip;               % [deg] fault dip
            strike = params.dip_azi - 90;   % [deg] fault strike    
            [sn0, tau0, ~] = project_stress_on_fault3D(SHdir, strike, dip, SH, Sh, Sv); % [MPa] total stresses
            self.tau0 = -tau0;
            self.sn0 = sn0;
        end
    end

end