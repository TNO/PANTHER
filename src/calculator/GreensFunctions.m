classdef GreensFunctions
    % GreensFunctions Returns Green's functions for a single HW or FW
    % compartment. (lenght(y), 1)
    properties
        Gxx_FW
        Gyy_FW
        Gxy_FW
        Gxx_HW
        Gyy_HW
        Gxy_HW
        Gnorm_FW
        Gshear_FW
        Gnorm_HW
        Gshear_HW
    end

    methods
        function self = GreensFunctions(yeval)
            % GreensFunctions computes Green's functions xx, yy, xy for FW
            % or HW compartments
            props = properties(self);
            for i = 1 : length(props)
                self.(props{i}) = zeros(length(yeval), 1);
            end
        end

        function self = green_FW(self, xeval, yeval, dip, h, t, w_FW, xcf, ycf)
            % green_FW computes Green's function for footwall compartment 
            % xeval     [array] evalution points in x
            % yeval     [array] evaluation points in y
            % dip       [deg] fault dip
            % h         [m] compartment or vertical slice height
            % t         [m] compartment or vertical slice throw 
            % w_FW      [m] width of footwall compartment (from xcf)
            % ycf       [m] y-coordinate mid depth of compartments
            % xcf       [m] x-coordinate at ycf
            [o, p, q, r, s] = self.get_corner_points_FW(dip, h, t, w_FW, xcf, ycf);
            for i = 1 : length(yeval)
                % xx components
                Green_triangle_xx = Gxx_triangle(o, p, r, s, xeval(i), yeval(i));
                Green_rectangle_xx = Gxx_rectangle(p, q, r, s, xeval(i), yeval(i));
                self.Gxx_FW(i) = Green_triangle_xx + Green_rectangle_xx;
                % yy components
                Green_triangle_yy = Gyy_triangle(o, p, r, s, xeval(i), yeval(i));
                Green_rectangle_yy = Gyy_rectangle(p, q, r, s, xeval(i), yeval(i));
                self.Gyy_FW(i) = Green_triangle_yy + Green_rectangle_yy;
                % xy components
                Green_triangle_xy = Gxy_triangle(o, p, r, s, xeval(i), yeval(i));
                Green_rectangle_xy = Gxy_rectangle(p, q, r, s, xeval(i), yeval(i));
                self.Gxy_FW(i) = Green_triangle_xy + Green_rectangle_xy;
            end
            [self.Gnorm_FW, self.Gshear_FW] = self.transform_principal_to_fault(self.Gxx_FW, self.Gyy_FW, self.Gxy_FW, dip);
            [self.Gnorm_FW] = self.subtract_2pi_FW(xeval, yeval, dip, h, t, ycf);% and add a criterium to check if the FW dP and dT are not 0;    
        end

        function self = green_HW(self, xeval, yeval, dip, h, t, w_HW, xcf, ycf)
            % green_HW computes Green's function for hanginwall compartment 
            % xeval     [array] evalution points in x
            % yeval     [array] evaluation points in y
            % dip       [deg] fault dip
            % h         [m] compartment or vertical slice height
            % t         [m] compartment or vertical slice throw 
            % w_HW      [m] width of hangingwall compartment (from xcf)
            % ycf       [m] y-coordinate mid depth of compartments
            % xcf       [m] x-coordinate at ycf
            [o, p, q, r, s] = self.get_corner_points_HW(dip, h, t, w_HW, xcf, ycf);
            for i = 1 : length(yeval)
                % NB the HW triangle is rotated. s<->r. 
                % xx components
                Green_triangle_xx = Gxx_triangle(o, q, s, r, xeval(i), yeval(i));
                Green_rectangle_xx = Gxx_rectangle(p, q, r, s, xeval(i), yeval(i));
                self.Gxx_HW(i) = Green_triangle_xx + Green_rectangle_xx;
                % yy components
                Green_triangle_yy = Gyy_triangle(o, q, s, r, xeval(i), yeval(i));
                Green_rectangle_yy = Gyy_rectangle(p, q, r, s, xeval(i), yeval(i));
                self.Gyy_HW(i) = Green_triangle_yy + Green_rectangle_yy;
                % xy components
                Green_triangle_xy = Gxy_triangle(o, q, s ,r, xeval(i), yeval(i));
                Green_rectangle_xy = Gxy_rectangle(p, q, r, s, xeval(i), yeval(i));
                self.Gxy_HW(i) = Green_triangle_xy + Green_rectangle_xy;
            end
            [self.Gnorm_HW, self.Gshear_HW] = self.transform_principal_to_fault(self.Gxx_HW, self.Gyy_HW, self.Gxy_HW, dip);
         %   [self.Gnorm_HW] = self.subtract_2pi_HW(xeval, yeval, dip, h, t, ycf);
        end

        function [o, p, q, r, s] = get_corner_points_FW(~, dip, h, t, w_FW, xcf, ycf)
            b = (h + t)/2;
            a = (h - t)/2;
            x_offset = 1/tan(dip*pi/180);
            o = xcf - a*x_offset;
            p = xcf + (b*x_offset);
            q = xcf + w_FW + (0.5*t*x_offset);
            r = ycf - a;
            s = ycf + b;
        end

        function [o, p, q, r, s] = get_corner_points_HW(~, dip, h, t, w_HW, xcf, ycf)
            b = (h + t)/2;
            a = (h - t)/2;
            x_offset = 1/tan(dip*pi/180);
            o = xcf + a*x_offset;
            p = xcf - w_HW - (0.5*t*x_offset);  % left corner
            q = xcf - (b*x_offset);             % right corner
            r = ycf - b;                        % top y
            s = ycf + a;                        % bottom y
        end

        function [normal, shear] = transform_principal_to_fault(~, xx, yy, xy, dip)
            normal = xx.*(sin(dip*pi/180)).^2 + yy .* (cos(dip*pi/180)).^2 - 2.*xy.*sin(dip*pi/180).*cos(dip*pi/180);
            shear = xy.*((sin(dip*pi/180)).^2 - (cos(dip*pi/180)).^2) + (xx - yy).*sin(dip*pi/180).*cos(dip*pi/180);
        end

        function [Gnorm_FW] = subtract_2pi_FW(self, xeval, yeval, dip, h, t, ycf)
            xfault = yeval/(tan(dip*pi/180));
            dx = xeval - xfault;
            b = (h + t)/2 ;     % with respect to ycf (= mid depth of compartment)
            a = (h - t)/2 ;     % with respect to ycf (= mid depth of compartment)   
            Gnorm_FW = self.Gnorm_FW;
            for i = 1 : length(dx)
               if dx(i) >= 0 && and((yeval(i) - ycf) <= b, (yeval(i) - ycf) >= -a)
                   Gnorm_FW(i) = Gnorm_FW(i) - 2 *pi;
               end
            end
        end

        function [Gnorm_HW] = subtract_2pi_HW(self, xeval, yeval, dip, h, t, ycf)
            xfault = yeval/(tan(dip*pi/180));
            dx = xeval - xfault;
            b = (h + t)/2 + ycf;
            a = (h - t)/2 + ycf;
            Gnorm_HW = self.Gnorm_HW;
            for i = 1 : length(dx)
               if dx(i) <= 0 && and(yeval(i) <= a, yeval(i) >= -b)
                   Gnorm_HW(i) = Gnorm_HW(i) - 2 *pi;
               end
            end
        end

    end

end