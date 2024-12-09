function [y_top] = get_top_y(thickness, throw, side)
    % function to compute the base top of the FW or HW compartment, with
    % resepct to the center depth (depth_ymid)
    % thick = [double] reservoir thickness
    % throw = [double] reservoir throw
    % side = [string] side, FW: footwall and HW: hanging wall
    if strcmp(p_side, 'FW')
        y_top = (thick + throw)/2;
    elseif strcmp(p_side, 'HW')
        y_top = (thick - throw)/2; 
    else
        error('Incorrect side indicator entered, should be FW or HW');
    end
end