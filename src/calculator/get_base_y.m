function [y_base] = get_base_y(thick, throw, side)
    % function to compute the base depth of the FW or HW compartment, with
    % resepct to the center depth (depth_ymid)
    % thick = [double] reservoir thickness
    % throw = [double] reservoir throw
    % side = [string] side, FW: footwall and HW: hanging wall
    if strcmp(side, 'FW')
        y_base = -(thick - throw)/2;
    elseif strcmp(side, 'HW')
        y_base = -(thick + throw)/2;
    else
        error('Incorrect side indicator entered, should be FW or HW');
    end
end