% set the vertical spacing along the fault
% TODO consider adjusting for fault dip to assure constant resolution along
% the fault 

function [y] = get_fault_y(dy, y_extent)
    
    ny = 1 + 2*floor(y_extent/dy);
    dy_2 = 2*y_extent/ny;           % Rescale dy
    y = zeros(ny, 1);
    
    for i = 1 : ny
        y(i) = y_extent - (i-1)*dy_2;
    end

end