% Green's functions for stress in y (vertical). 
% The reservoir is modeled through a combination of a rectangle and triangle 
% Can be for reservoirs on two sides of the fault or one side
% c = 0 --> right side only
% d = 0 --> left side only

function G_yy = Gyy_combined(a,b,c,d,theta,x,y)

num_offset = 1e-9; % small number to add to x and y to avoid singularities and division by 0) [-]  
% why time abs(d)?
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
    Gyy_left_rectangle = Gyy_rectangle(-c,-f,-b,a,x,y);
    Gyy_left_triangle  = Gyy_triangle(-e,f,-a,b,-x,-y);
    Gyy_right_triangle = Gyy_triangle(-e,f,-a,b,x,y);
    Gyy_right_rectangle = Gyy_rectangle(f, d,-a,b,x,y);
elseif c==0 && ~d==0
    % right side only
    Gyy_left_rectangle = 0;
    Gyy_left_triangle  = 0;
    Gyy_right_triangle = Gyy_triangle(-e,f,-a,b,x,y);
    Gyy_right_rectangle = Gyy_rectangle(f, d,-a,b,x,y);
elseif ~c==0 && d==0
     % left side only
    Gyy_left_rectangle = Gyy_rectangle(-c,-f,-b,a,x,y);
    Gyy_left_triangle  = Gyy_triangle(-e,f,-a,b,-x,-y);
    Gyy_right_triangle = 0;
    Gyy_right_rectangle = 0;
end

G_yy = Gyy_left_rectangle + Gyy_left_triangle + Gyy_right_rectangle + Gyy_right_triangle;

end


