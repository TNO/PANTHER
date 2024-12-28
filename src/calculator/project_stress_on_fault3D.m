function [s_normal, tau_dip, tau_strike] = project_stress_on_fault3D(SHdir, strike, dip, SHall, Shall, Svall)
    
    % take Sh SH and SV as the principal stress    
    % SHdir, strike, and dip can be the a constant, different array lengths of
    % SH, Sv, and Sh allowed
    
    
    s_normal = zeros(size(Shall));      % total normal stress
    tau_dip = zeros(size(Shall));       % shear stress along dip
    tau_strike = zeros(size(Shall));    % shear stress along strike

    for i = 1 : length(Shall)
        Sh = Shall(i);
        SH = SHall(i);
        Sv = Svall(i);
        if length(SHdir) == length(Shall)
            SHdir_r = SHdir(i)*pi/180;
        else
            SHdir_r = SHdir(1)*pi/180;
        end

        if length(strike) == length(Shall)
            dstrike_r = SHdir_r - strike(i)*pi/180; % angle between strike and SHmax direction
        else
            dstrike_r = SHdir_r - strike(1)*pi/180;
        end

        if length(dip) == length(Shall)
            dip_r = dip(i)*pi/180;
        else
            dip_r = dip(1)*pi/180;
        end

        if Sh > SH
            error("Sh is larger than SH");
        end
    
        % stress tensor of principal stresses
        Scomp = sort([Sh, SH, Sv], 'descend');
        SP = zeros(3,3);
        SP(1,1) = Scomp(1); % s1
        SP(2,2) = Scomp(2); % s2
        SP(3,3) = Scomp(3); % s3
       
       
        % unit vectors defining fault plane
        % n_n unit vector normal to fault plane 
        % normal faulting
         if ((Sh <= SH) && (Sh < Sv))
            n_n(1) = -cos(dip_r);                   % angle normal vector with sv (here:s1)
            n_n(2) = -sin(dstrike_r)*sin(dip_r);    % angle normal vector with sH (here: s2)
            n_n(3) = cos(dstrike_r)*sin(dip_r);     % angle normal vector with sh (here: s3)
            n_s(1) = 0;                             % angle strike vector with sv  (90deg)
            n_s(2) = cos(dstrike_r);                % angle strike vector with sH 
            n_s(3) = sin(dstrike_r);                % angle strike vector with sh
            n_d(1) = sin(dip_r);                    % angle dip vector with sv
            n_d(2) = -sin(dstrike_r)*cos(dip_r);    % angle dip vector with sH
            n_d(3) = cos(dstrike_r)*cos(dip_r);     % angle dip vector with sh
         % strike-slip faulting, Sv = s2
        elseif (Sh <= SH) && (SH > Sv)
            n_n(1) = -sin(dstrike_r)*sin(dip_r);    % angle normal vector with sH (here:s1)
            n_n(2) = cos(dstrike_r)*sin(dip_r);     % angle normal vector with sv (here:s2)
            n_n(3) = -cos(dip_r);                   % angle normal vector with sh (here:s3)
            % unit vector along strike
            n_s(1) = cos(dstrike_r);
            n_s(2) = sin(dstrike_r);
            n_s(3) = 0;
            % unit vector along dip
            n_d(1) = -sin(dstrike_r)*cos(dip_r);
            n_d(2) = cos(dstrike_r)*cos(dip_r);
            n_d(3) = sin(dip_r);
        % thrust faulting
        elseif (Sh <= SH) && (Sh > Sv)
            n_n(1) = -sin(dstrike_r)*sin(dip_r);
            n_n(3) = cos(dstrike_r)*sin(dip_r);
            n_n(2) = -cos(dip_r);
            % unit vector along strike
            n_s(1) = cos(dstrike_r);
            n_s(3) = sin(dstrike_r);
            n_s(2) = 0;
            % unit vector along dip
            n_d(1) = -sin(dstrike_r)*cos(dip_r);
            n_d(3) = cos(dstrike_r)*cos(dip_r);
            n_d(2) = sin(dip_r);
         end
        
        % fault stresses
        tF = SP*n_n';       % fault traction vector
        s_normal(i) = dot(tF,n_n);   % normal stress
        tau_dip(i) = dot(tF,n_d);    % shear stress in dip direction
        tau_strike(i) = dot(tF,n_s); % shear stress along strike
    
        % tau_max(i) = (norm(tF)^2 - norm(s_normal(i))^2)^0.5;
        % tau_max(i) = (tau_dip(i)^2 + tau_strike(i)^2)^0.5;

        if abs(tau_dip(i)) < 1e-12
            tau_dip(i) = 0;
        end
    
    end
    % https://dnicolasespinoza.github.io/node38.html#table:RPGsummary

end

