function [next_to_FW, next_to_HW, next_to_res] = is_adjacent_to_reservoir(y, thickness, offset)
    % is_adjacent_to_reservoir Identifies which fault depth intervals 
    % are adjacent to one or both of the reservoir compartments
    % INPUTS
    % y             = depth relative to depth_mid
    % thickness     = reservoir thickness
    % offset        =  vertical reservoir offset
    % RETURNS
    % next_to_HW    = logical array length(y) with 1 indicating adjacent to HW
    % next_to_FW    = logical array length(y) with 1 indicating adjacent to FW
    % next_to_res    = logical array length(y) with 1 indicating adjacent
    % to FW or HW

    a = (thickness - offset)/2;  % +ve top left compartment, -ve base right
    b = (thickness + offset)/2;  % 
    top_FW = b; 
    base_FW = -a;
    top_HW = a; 
    base_HW = -b;
    
    next_to_HW = (y >= base_HW & y <= top_HW);
    next_to_FW = (y >= base_FW & y <= top_FW);
    next_to_res = or(next_to_HW, next_to_FW);

end