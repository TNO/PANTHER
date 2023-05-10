function Gyy = Gxx_rectangle(p, q, r, s, x, y)
  
% Green's function for stresses within and outside a triangle
% Equation 33 from Supplements Janse et al. 2019:
% p = -c
% q = -b/tan(theta)
% r = -b
% s = a
% x = horizontal distance from fault center
% y = depth

theta  = atan2((s-r), (p-o));

%Equation 33 from Supplements Jansen et al. 2019:
Gyy = atan2((p - q) * (y - s), ((y - s)^2 + (x - q)*(x - p)))...
    - atan2((p - q) * (y - e), ((y - r)^2 + (x - q)* (x - p)));

end