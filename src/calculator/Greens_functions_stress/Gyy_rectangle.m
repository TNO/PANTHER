function Gyy = Gyy_rectangle(p, q, r, s, x, y)
% TODO allow inf width

% Green's function for stresses within and outside a rectangle, y
% Equation 33 from Supplements Jansen et al. 2019:

% theta = dip angle [rad]
% p and q are the corners of the rectangle in x
% r and s are the corner of the rectangle in y
% p = left edge =   -c (LEFT block)   or        b/tan(theta) (RIGHT)    
% q = right edge =  -b/tan(theta) (LEFT)        d (RIGHT block)
% r = bottom edge = -b (LEFT)                  -a (RIGHT)
% s = top edge =     a (LEFT)                   b (RIGHT)
% x = horizontal distance from fault center
% y = depth

% quick fix for infinite. TODO derive exact equation for rectangle
if q == inf
    q = 1e15;
end
if q == -inf
    q = -1e15;
end
if p == inf
    p = 1e15;
end
if p == -inf
    p = -1e15;
end


%Equation 33 from Supplements Jansen et al. 2019:
Gyy = atan2((p - q) * (y - s), ((y - s)^2 + (x - q)*(x - p)))...
    - atan2((p - q) * (y - r), ((y - r)^2 + (x - q)* (x - p)));

end