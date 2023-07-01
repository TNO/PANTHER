function [T_reservoir] = calc_dT_diffusion(y, y_top, y_base, time_steps, T_reservoir, diffusivity)
    % calc_dT_diffusion computes temperature diffusion to seal and base, based
    % on the reservoir temperature and the distance to reservoir top or base
    % uniform temperature inside the reservoir. steady-state diffusion in
    % assumed at each t, into both the seal and underburden, which have the
    % same diffusivity. Diffusion is purely vertical. Mossop (2001). 
    % INPUT
    % y             [m] depth, relative to mid reservoir depth 
    % time          [yr] operation time
    % T_reservoir   [MPa] reservoir temperature (size length(y) x length(t))
    % i_reservoir   [-] indices of the reservoir section
    % OUTPUT
    % T_reservoir [MPa] T_reservoir array updated with diffusion
    % temperatures
    
    % top indices reservoir section
    top_i = find(y < y_top, 1, 'first');  % index of first non-zero pressure value (top reservoir )
    base_i = find(y_base < y , 1, 'last');  % index of last non-zero pressure value (base reservoir)
    
    time_s_all = time_steps * 365.25* 24 * 60 *60;
    % diffusion into the seal
    diff_y = y - y_top; 
    T_temp = T_reservoir(top_i, :) .* erfc(diff_y./(2*sqrt(diffusivity*time_s_all')));
    T_reservoir(1:top_i-1, :) = T_temp(1:top_i-1, :);
    % diffusion into the underburden
    diff_y = y_base - y; 
    T_temp = T_reservoir(base_i, :) .* erfc(diff_y./(2*sqrt(diffusivity*time_s_all')));
    T_reservoir(base_i+1:end, :) = T_temp(base_i+1:end, :);
    
end