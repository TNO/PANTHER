classdef InitialStress
    % InitialPressure sets the initial pore pressure in the reservoir
    % compartments and in the fault. For different density reservoir fluids
    % (e.g. methane) this density is used to compute the pressure (i.e. the gas pressure gradient) 

    properties
        sn0 double       % [MPa] initial total normal stress
        tau0 double    % [MPa] initial shear stress
    end

    methods
        function self = InitialStress(y, params)
            % Initializes initial stress
            % TODO allow for non-uniform initial stress
            y_abs = y + params.depth_mid;
            Sv = -(y_abs/1000) .* params.sv_grad + params.sv_offset;    % [MPa] vertical stress
            Sh = params.shsv .* Sv;         % [MPa] total minimum horizontal stress
            SH = Sh.*params.sHsh;           % [MPa] total maximum horizontal stress
            SHdir = params.sH_dir;          % [deg] direction of maximum horizontal stress
            dip = params.dip;               % [deg] fault dip
            strike = params.dip_azi - 90;   % [    
            [sn0, tau0, ~] = project_stress_on_fault3D(SHdir, strike, dip, SH, Sh, Sv); % [MPa] total stresses
            self.tau0 = -tau0;
            self.sn0 = sn0;
        end
    end

end