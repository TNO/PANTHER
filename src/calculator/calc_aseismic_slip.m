function [tau_slip, islip, slip] = calc_aseismic_slip(y, tau, strength, mu_II, K)
    
    % find indices of slipping cells
    islip = find(tau >= strength);      % find indices of slipping cells
    tau_slip = tau; 
% 
%     % influence function, from Pablo's code td_solve 
%     Kline = -nx/(2*pi*(y(1)-y(end))) ./ ( [0:nx-1]'.^2-0.25 );
%     Ko = toeplitz(Kline);                   
%     K = mu_II*Ko;
%     K(and(K<0,K>(min(min(K)/10000)))) = 0;    % set very small changes to 0 to avoid continued interactions of the two peaks to 
%     
    slip = zeros(size(tau));
    while 1
        % identify number of zones where tau>strength
        tauexc = -(tau(islip) - strength(islip));   % excess shear stress in slip zones
        [tau_red, slip_red] = td_solve(islip, tauexc, 0, K);  % tauexc: stress drop for islip, stress increase outside
        tau = tau + tau_red;                         % update shear stress with redistributed stress
        slip = slip + slip_red;
        islipred = find(tau >= strength);           % find slipping cells for redistributed stress state
        if isequal(islipred,islip)
            break      % if no new cells reached the strength, stop redistribution
        else
             islip = islipred;                    % update cell indices of slipping cells
        end
    end
    tau_slip = tau;

end