function Gyy = Gyy_triangle(o, p, r, s, x, y)
% Green's function for stresses within and outside a triangle
% Equation 36 from Supplements Jansen et al. 2019:
% theta = dip angle [rad]
% G_yy_2  = g_yy_int_triang(-e,f,-a,b,-x,-y);
% G_yy_3 = g_yy_int_triang(-e,f,-a,b,x,y);
% -e,f,-a,b,
% same input for the triangle, except x and y sign (180 degree turn)
% o = -a/tan(theta)     x coor of angle at the fault zone (topmost for t>0)
% p = b/tan(theta)      x coor of right angle 
% r = -a                y coor of right angle
% s = b                 y coor of angle at the fault zone
% x = horizontal distance from fault center
% y = depth

    theta  = atan2((s-r), (p-o));

%     %Equation 36 from Supplements Jansen et al. 2019:
%     Gyy    =    (1/4)*sin(2*theta)*log(((y-p)^2 + (x - p*tan(theta))^2)...
%         /((x - o)^2 + (y - o*tan(theta))^2))...
%         -atan2((o-p)*(y-r),((y-r)^2 + (x - p)*(x - o)))...
%         -(cos(theta))^2*atan2((p-o)*(y-x*tan(theta)),(x^2+y^2+o*p*(sec(theta))^2-(p+o)*(x+y*tan(theta))));

    Gyy = 0.5 * log(((x - p)^2 + (y - p*tan(theta))^2) ...
        / ((x - o)^2 + (y + o*tan(theta))^2)) * sin(theta)*cos(theta) ...
        - atan2((o - p)*(y - r) ,((y - r)^2 + (x - p)*(x - o)))...
        - atan2((p - o)*(y - x*tan(theta)), ...
        (x^2 + y^2 + o*p*(sec(theta))^2 - (p + o)*(x + y*tan(theta)))) * (cos(theta))^2;

end