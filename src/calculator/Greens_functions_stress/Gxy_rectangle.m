function Gxy = Gxy_rectangle(p, q, r, s, x, y)
  
% Green's function for stresses within and outside a rectangle, shear
% Equation 34 from Supplements Jansen et al. 2019:

% theta = dip angle [rad]
% p and q are the corners of the rectangle in x
% r and s are the corner of the rectangle in y
% p = left edge =   -c (LEFT block)   or        b/tan(theta) (RIGHT)    
% q = right edge =  -b/tan(theta) (LEFT)        d (RIGHT block)
% r = bottom edge = -b (LEFT)                  -a (RIGHT)
% s = top edge =     a (LEFT)                   b (RIGHT)
% x = horizontal distance from fault center
% y = depth

% theta  = atan2((s-r), (p-o));

% Equation 34 from Supplements Jansen et al. 2019:
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


Gxy = 0.5 * log( ((x - q)^2 + (y - s)^2) * ((x - p)^2 + (y - r)^2) ...
        /  (((x - q)^2 + (y - r)^2) * ((x - p)^2 + (y - s)^2)) );
    
end