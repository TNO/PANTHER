function Gxy = Gxy_triangle(o, p, r, s, x, y)
% Green's function for stresses within and outside a triangle
% Equation 37 from Supplements Jansen et al. 2019:
% theta = dip angle [rad]
% same input for the triangle, except x and y sign (180 degree turn)
% o = -a/tan(theta)     x coor of angle at the fault zone (topmost for t>0)
% p = b/tan(theta)      x coor of right angle 
% r = -a                y coor of right angle
% s = b                 y coor of angle at the fault zone
% x = horizontal distance from fault center
% y = depth
% G_xy_3 = g_xy_int_triang(-e,f,-a,b,x,y);

    theta  = atan2((s-r), (p-o));   % dip angle (also allowing for theta>0.5pi

    %Equation 37 from Supplements Jansen et al. 2019:
    Gxy = 0.5 * log(((x - p)^2 +(y - s)^2)/((x - p)^2 + (y - r)^2)) ...
        + 0.5 * (sin(theta))^2 * log(((y - r)^2 + (x - r*cot(theta))^2) ...
        / ((y - s)^2  + (x - s * cot(theta))^2)) ...
    + atan2(((s - r)*(y*cot(theta) - x)), ...
        (x^2 + y^2 + r*s*(csc(theta))^2 - (r + s)*(y + x*cot(theta)))) * sin(theta) * cos(theta);     
    
% Gxy    =    (1/2)*log(((x-p)^2+(y-s)^2)/((x-p)^2+(y-r)^2))...
%     +(1/2)*(sin(theta))^2*log(((y-r)^2+(x-r*cot(theta))^2)...
%     /((y-s)^2+(x-s*cot(theta))^2))+sin(theta)*cos(theta)...
%     *atan2(((s-r)*(y*cot(theta)-x)),(x^2+y^2+r*s*(csc(theta))^2-(r+s)*(y+x*cot(theta))));
%     
%     + (x - r*cot(theta))^2) ...
%         / ((y - s)^2 + (x - s*cot(theta))^2)) * sin(theta)*cos(theta) ...
%         + atan2((r - s)*(x - p) ,((x - p)^2 + (y - r)*(y - s)))...
%         - atan2((r - s)*(y*cot(theta) - x), ...
%         (x^2 + y^2 + r*s*(csc(theta))^2 - (r + s)*(y + x*cot(theta)))) * (sin(theta))^2;

end