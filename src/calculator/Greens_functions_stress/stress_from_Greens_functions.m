function [G_norm, dsigma_n, dtau, x, y] = stress_from_Greens_functions(varin)


a = varin.a;        % y distance to inner reservoir corners, from center
b = varin.b;        % y distance to outer reservoir corners, from center
c = varin.c;        % width of left side
d = varin.d;        % width of right side
C = varin.C;        % coefficient for displacements
g = varin.g;        % gravitational constant
ny = varin.ny;
theta = varin.theta;
y_max = varin.y_max;
y_min = varin.y_min;
dx = varin.dx;

small = 1e-9;

% Computes stresses at the fault or at a line parallel to the fault

% initialize arrays for Gnorm and stress changes
Dy_2 = (y_max - y_min)/ny;
y = zeros(ny,1);
G_norm_minus2pi = zeros(ny,1);
G_norm = zeros(ny,1);
G_shear = zeros(ny,1);
dtau = zeros(ny,1); 
dsigma_n = zeros(ny,1);

for i = 1 : ny
        y(i) = y_min + (i-1)*Dy_2;
        x = y(i)*cot(theta) + small + dx;
        % Gxx, Gyy, Gxy 
        G_xx = Gxx_combined(a,b,c,d,theta,x,y(i));
        G_yy = Gyy_combined(a,b,c,d,theta,x,y(i));
        G_xy = Gxy_combined(a,b,c,d,theta,x,y(i));
        % transform x y to fault-perpendicular and fault-parallel
        G_norm(i)  = G_xx * (sin(theta))^2 + G_yy * (cos(theta))^2 - 2*G_xy * sin(theta) * cos(theta);
        % subtract factor 2 pi from solution, only for right side of the
        % fault, because we are always considering x > 0 
        % see Jansen supplements text after Eq 51
        if  (y(i) >= -a && y(i) <= b && d~=0)      % rigth side only
            G_norm_minus2pi(i) = G_norm(i)- 2*pi;
        elseif (y(i) >= - b && y(i) <= a && d==0) % left side only 
            G_norm_minus2pi(i) = G_norm(i);    
        else
            G_norm_minus2pi(i) = G_norm(i);
        end
        G_shear(i) = G_xy * ( (sin(theta))^2 - (cos(theta))^2 ) + (G_xx - G_yy) * sin(theta) * cos(theta);
end
dsigma_n = C * G_norm_minus2pi;
dtau = C * G_shear;

end 