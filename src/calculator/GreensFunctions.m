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
        y_correction = 0    % (length(y) x 1) correction applied to points coinciding with a reservoir boundary
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
            % use the uncorrected y to identify whether a calculation point
            % lies within the reservoir or at the reservoir boundary
            yeval_uncorrected = yeval - self.y_correction;
            for i = 1 : length(dx)
               if dx(i) >= 0 && and((yeval_uncorrected(i) - ycf) <= b, (yeval_uncorrected(i) - ycf) >= -a)
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
            % use the uncorrected y to identify whether a calculation point
            % lies within the reservoir or at the reservoir boundary
            yeval_uncorrected = yeval - self.y_correction;
            for i = 1 : length(dx)
               if dx(i) <= 0 && and(yeval_uncorrected(i) <= a, yeval_uncorrected(i) >= -b)
                   Gnorm_HW(i) = Gnorm_HW(i) - 2 *pi;
               end
            end
        end

    end

    methods (Static)

        function GF = initialize(params, y, dx, variable_PT, variable_dip)
            % initialize  Build a cell array of GreensFunctions objects for
            % a given fault geometry, choosing the appropriate mode:
            %   uniform    — single GF shifted along-depth  (fastest)
            %   variable PT — one GF convolved by shift vector (intermediate)
            %   variable dip — separate GF at every depth cell (slowest)
            %
            % INPUT
            % params        PantherMember with fault/reservoir geometry
            % y             [m] depth array w.r.t. y_mid
            % dx            [m] distance from fault in x
            % variable_PT   logical — true when P or T vary with depth
            % (vertical slice can be translated)
            % variable_dip  logical — true when dip varies with depth
            % (separate Greens's function at each depth)
            %
            % OUTPUT
            % GF   cell(1,1) for uniform case; cell(length(y),1) otherwise

            correction_value = 1e-3;
            if dx == 0 || dx == params.width_FW || dx == -params.width_HW
                dx = dx + correction_value;
            end
            reservoir_boundaries = [params.y_FW_top, params.y_HW_top, params.y_HW_base, params.y_FW_base];
            y_correction = zeros(size(y));
            for i = 1 : length(reservoir_boundaries)
                if ismember(y, reservoir_boundaries(i))
                    if i == 1 || i == 2
                        y_correction(ismember(y, reservoir_boundaries(i))) = -correction_value;
                    else
                        y_correction(ismember(y, reservoir_boundaries(i))) = correction_value;
                    end
                end
            end
            y = y + y_correction;
            xeval = y ./ (tan(params.dip * pi / 180)) + dx;

            if ~variable_PT && ~variable_dip
                % --- uniform case ---
                GF{1} = GreensFunctions(y);
                GF{1}.y_correction = y_correction;
                if params.width_FW > 0
                    GF{1} = GF{1}.green_FW(xeval, y, params.dip, params.thick, params.throw, params.width_FW, 0, 0);
                end
                if params.width_HW > 0
                    GF{1} = GF{1}.green_HW(xeval, y, params.dip, params.thick, params.throw, params.width_HW, 0, 0);
                end

            elseif variable_PT && ~variable_dip
                % --- variable P/T, uniform dip: shift a single GF ---
                slice_thick = y(1) - y(2);
                y2     = [y; y(1:end-1) + (y(end) - y(1)) - slice_thick];
                xeval2 = y2 / (tan(params.dip * pi / 180)) + dx;
                i_mid  = ceil(length(y2) / 2);
                slice_y = y2(i_mid);
                slice_x = slice_y / (tan(params.dip * pi / 180));
                slice_throw = 0;
                greens_f = GreensFunctions(y2);
                if params.width_FW > 0
                    greens_f = greens_f.green_FW(xeval2, y2, params.dip, slice_thick, slice_throw, params.width_FW, slice_x, slice_y);
                end
                if params.width_HW > 0
                    greens_f = greens_f.green_HW(xeval2, y2, params.dip, slice_thick, slice_throw, params.width_HW, slice_x, slice_y);
                end
                GF = cell(length(y), 1);
                for j = 1 : length(y)
                    GF{j} = greens_f;
                    GF{j}.Gnorm_FW  = GF{j}.Gnorm_FW(i_mid-j+1 : 2*i_mid-j);
                    GF{j}.Gnorm_HW  = GF{j}.Gnorm_HW(i_mid-j+1 : 2*i_mid-j);
                    GF{j}.Gshear_FW = GF{j}.Gshear_FW(i_mid-j+1 : 2*i_mid-j);
                    GF{j}.Gshear_HW = GF{j}.Gshear_HW(i_mid-j+1 : 2*i_mid-j);
                end

            else
                % --- variable dip: separate GF per depth cell (slowest) ---
                slice_thick = y(1) - y(2);
                slice_throw = 0;
                GF = cell(length(y), 1);
                for j = 1 : length(y)
                    if variable_dip
                        dip_j = params.dip(j);
                    else
                        dip_j = params.dip;
                    end
                    slice_y = y(j);
                    slice_x = slice_y / (tan(dip_j * pi / 180));
                    GF{j} = GreensFunctions(y);
                    GF{j} = GF{j}.green_FW(xeval, y, dip_j, slice_thick, slice_throw, params.width_FW, slice_x, slice_y);
                    GF{j} = GF{j}.green_HW(xeval, y, dip_j, slice_thick, slice_throw, params.width_HW, slice_x, slice_y);
                end
            end
        end

    end

end