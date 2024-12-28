function Gxy = Gxy_combined(a,b,c,d,theta,x,y)
%
% Green's functions for shear stresses . 
% The reservoir is modeled through a combination of a rectangle and triangle 
% Can be for reservoirs on two sides of the fault or one side
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

% Compute integrals:
f = b/tan(theta);
e = a/tan(theta);

if and(~(c==0), ~d==0)
    Gxy_left_rectangle = Gxy_rectangle(-c,-f,-b,a,x,y);
    Gxy_left_triangle  = Gxy_triangle(-e,f,-a,b,-x,-y);
    Gxy_right_triangle = Gxy_triangle(-e,f,-a,b,x,y);
    Gxy_right_rectangle = Gxy_rectangle(f, d,-a,b,x,y);
elseif c==0 && ~d==0
    % right side only
    Gxy_left_rectangle = 0;
    Gxy_left_triangle  = 0;
    Gxy_right_triangle = Gxy_triangle(-e,f,-a,b,x,y);
    Gxy_right_rectangle = Gxy_rectangle(f, d,-a,b,x,y);
elseif ~c==0 && d==0
     % left side only
    Gxy_left_rectangle = Gxy_rectangle(-c,-f,-b,a,x,y);
    Gxy_left_triangle  = Gxy_triangle(-e,f,-a,b,-x,-y);
    Gxy_right_triangle = 0;
    Gxy_right_rectangle = 0;
end

Gxy = Gxy_left_rectangle + Gxy_left_triangle + Gxy_right_rectangle + Gxy_right_triangle;

end



