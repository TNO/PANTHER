function Gxx = Gxx_combined(a,b,c,d,theta,x,y)
    % Gxx_combined combines Green's functions for different reservoir compartments.
    %
    % This function combines Green's functions for different reservoir compartments.
    % The reservoir is modeled through a combination of a rectangles and triangles.
    %
    % Parameters:
    %   a - Vertical distance to inner reservoir corners, from center (y = 0)
    %   b - Vertical distance to outer reservoir corners, from center (y = 0)
    %   c - Length of the left side of the reservoir
    %   d - Length of the right side of the reservoir 
    %   theta - Angle of the fault [radians]
    %   x - x-coordinate of evaluation point(s)
    %   y - y-coordinate of evaluation point(s)
    %
    % Returns:
    %   Gxx - Combined Green's function value
    %
    % Notes:
    %   - If c = 0, the reservoir is on the right side only.
    %   - If d = 0, the reservoir is on the left side only.
    %   - A small number is added to x and y to avoid singularities and division by zero.
    %
    % Example:
    %   Gxx = Gxx_combined(10, 20, 5, 5, pi/4, 2, 3);
    
    % Combines Green's functions for the different reservoir compartments. 
    % The reservoir is modeled through a combination of a rectangle and triangle 
    % The function can be used for reservoirs on two sides of the fault or 
    % one side
    % c = 0 --> right side only
    % d = 0 --> left side only
    
    num_offset = 1e-9; % small number to add to x and y to avoid singularities and division by 0) [-]  
    if x == 0 || x == d
        x = x + num_offset*abs(d);
        x = x + num_offset;
    end
    if x == -c 
        x = x - num_offset*abs(c);
        x = x - num_offset;
    end
    if y == a || y == b
        y = y + num_offset*abs(b-a);
    end
    if y == -a || y == -b
        y = y - num_offset*abs(b-a);
    end
    
    % Compute integrals of individual components
    f = b/tan(theta);       % offset in x 
    e = a/tan(theta);
    
    if and(~(c==0), ~d==0)
        Gxx_left_rectangle = Gxx_rectangle(-c,-f,-b,a,x,y);
        Gxx_left_triangle  = Gxx_triangle(-e,f,-a,b,-x,-y);
        Gxx_right_triangle = Gxx_triangle(-e,f,-a,b,x,y);
        Gxx_right_rectangle = Gxx_rectangle(f, d,-a,b,x,y);
    elseif c==0 && ~d==0
        % right side only
        Gxx_left_rectangle = 0;
        Gxx_left_triangle  = 0;
        Gxx_right_triangle = Gxx_triangle(-e,f,-a,b,x,y);
        Gxx_right_rectangle = Gxx_rectangle(f, d,-a,b,x,y);
    elseif ~c==0 && d==0
         % left side only
        Gxx_left_rectangle = Gxx_rectangle(-c,-f,-b,a,x,y);
        Gxx_left_triangle  = Gxx_triangle(-e,f,-a,b,-x,-y);
        Gxx_right_triangle = 0;
        Gxx_right_rectangle = 0;
    end
    
    Gxx = Gxx_left_rectangle + Gxx_left_triangle + Gxx_right_rectangle + Gxx_right_triangle;

end





