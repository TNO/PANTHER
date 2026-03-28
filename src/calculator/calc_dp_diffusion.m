function [p_reservoir] = calc_dp_diffusion(y, y_top, y_base, time_steps, p_reservoir, diffusivity)
    % calc_dp_diffusion computes pressure diffusion to seal and base, based
    % on the reservoir pressure and the distance to reservoir top or base
    % uniform pressure inside the reservoir. steady-state diffusion in
    % assumed at each t, into both the seal and underburden, which have the
    % same diffusivity
    % INPUT
    % y             [m] depth, relative to mid reservoir depth 
    % time          [yr] operation time relative to the first row
    % p_reservoir   [MPa] reservoir pressure (size length(y) x length(t))
    % i_reservoir   [-] indices of the reservoir section
    % OUTPUT
    % p_reservoir [MPa] p_reservoir array updated with diffusion pressures
    
    % top_i: index of the last calculation point within the reservoir (upwards) or at the reservoir top boundary 
    top_i = find(y <= y_top, 1, 'first');       
    % base_i: index of the last calculation point within the reservoir (downwards) or at the reservoir base boundary 
    base_i = find(y_base <= y , 1, 'last');     
    
    time_s_all = (time_steps - time_steps(1)) * 365.25* 24 * 60 *60;    % time in s, relative to calculation start
    % diffusion into the seal
    diff_y = y - y_top;     % distance from the top of the reservoir compartment
    p_ini = p_reservoir;    % pressure in the reservoir compartment of interest
    if ~isempty(top_i)
        % pressure_differential: pressure difference between reservoir and the seal
        pressure_differential = p_reservoir(top_i, :) - p_reservoir(top_i-1, :);
        % p_temp = p_reservoir(top_i, :) .* erfc(diff_y./(2*sqrt(diffusivity*time_s_all')));
        p_temp = pressure_differential .* erfc(diff_y./(2*sqrt(diffusivity*time_s_all')));
        p_reservoir(1:top_i-1, :) = p_reservoir(1:top_i-1, :) + p_temp(1:top_i-1, :);
    else 
        
    end
    % diffusion into the underburden
    if ~isempty(base_i)
        diff_y = y_base - y; 
        pressure_differential = p_reservoir(base_i, :) - p_reservoir(base_i + 1, :);
        p_temp = pressure_differential .* erfc(diff_y./(2*sqrt(diffusivity*time_s_all')));
        p_reservoir(base_i+1:end, :) = p_reservoir(base_i+1:end, :) + p_temp(base_i+1:end, :);
    end
end

% function [pressure] =  calc_pressure_diffusion_single_side(x, diffusivity, time, p_boundary, p_background)    
%     pressure = p_background + p_boundary .* erfc(x./(2*sqrt(diffusivity*time_s_all')));    
% end