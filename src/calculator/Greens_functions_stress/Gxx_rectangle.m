function Gxx = Gxx_rectangle(p, q, r, s, x, y)
  
% Green's function for stresses within and outside a rectangle, x
% Equation 32 from Supplements Jansen et al. 2019:

% theta = dip angle [rad]
% p and q are the corners of the rectangle in x (right and left)
% r and s are the corner of the rectangle in y (top and bottom)
% p = left edge =   -c (LEFT block)   or        b/tan(theta) (RIGHT)    
% q = right edge =  -b/tan(theta) (LEFT)        d (RIGHT block)
% r = bottom edge = -b (LEFT)                  -a (RIGHT)
% s = top edge =     a (LEFT)                   b (RIGHT)
% x = horizontal distance from fault center, evaluation point
% y = depth, evaluation point

% theta  = atan2((s-r), (p-o));

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


%Equation 32 from Supplements Jansen et al. 2019:

Y1 = (r - s) * (x - q) ;
X1 = ((x - q)^2 + (y - r)*(y - s));
Y2 = (r - s) * (x - p);
X2 = ((x - p)^2 + (y - r)*(y - s));

Gxx = atan2((r - s) * (x - q), ((x - q)^2 + (y - r)*(y - s))) ...
    - atan2((r - s) * (x - p), ((x - p)^2 + (y - r)*(y - s)));
% 
% Gxx = atan( Y1 / X1) + (pi/2)*sign(Y1)*(1-sign(X1)) ...
%     - atan( Y2 / X2) + (pi/2)*sign(Y2)*(1-sign(X2));

% strange behavior? Gxx increasing for larger values of p, but then beyond
% some value decreases again to zero before suddenly switching to pi/4
% atan2((r - s) * (x - p), ((x - p)^2 + (y - r)*(y - s)))
% p = -100, 0.0706, p = -1000, -0.1663, p = 1e-4, -0.024

% seems to give some error
% G_xx = atan2( (r - s)*(x - q) , ((x - q)^2 + (y - r)*(y - s)) ) ...
%      - atan2( (r - s)*(x - p) , ((x - p)^2 + (y - r)*(y - s)) );

end