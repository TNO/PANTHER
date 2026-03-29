classdef (HandleCompatible) ModelGeometry
    % ModelGeometry Class to represent the geometry of a 2D geological 
    % model consisting of two reservoir compartments offset by a fault.
    % This class contains properties and methods to calculate various
    % geometric parameters of the model, such as depths and indices of
    % hanging wall and footwall compartments.
    %
    % Properties:
    %   depth_mid   Mid depth of the model (must be negative)
    %   thick       Thickness of the model (must be positive)
    %   dip         [degrees] Dip angle of the model
    %   dip_azi     [degrees] Dip azimuth of the model
    %   throw       Throw of the model
    %   width_HW    Width of the hanging wall (default is Inf)
    %   width_FW    Width of the footwall (default is Inf)
    %
    % Methods:
    %   ModelGeometry - Constructor to initialize the model geometry
    %   y_HW_top - Returns depth of the top hanging wall, relative to mid depth
    %   y_FW_top - Returns depth of the top footwall, relative to mid depth
    %   y_HW_base - Returns depth of the base hanging wall, relative to mid depth
    %   y_FW_base - Returns depth of the base footwall, relative to mid depth
    %   i_HW_top - Returns index of the first element in the top hanging wall
    %   i_HW_base - Returns index of the last element in the base hanging wall
    %   i_FW_top - Returns index of the first element in the top footwall
    %   i_FW_base - Returns index of the last element in the base footwall
    %   get_y_base - Computes the base depth of the FW or HW compartment
    %   get_y_top - Computes the top depth of the FW or HW compartment
    %   i_FW - Returns indices of all elements in the footwall
    %   i_HW - Returns indices of all elements in the hanging wall
    %   i_reservoir - Returns indices of all elements in FW or HW

    properties
        depth_mid (1,1) double {mustBeNegative} = -3000
        thick (1,1) double {mustBePositive} = 150
        dip (:,1) double {mustBeInRange(dip, 0, 90)} = 60
        dip_azi (:,1) double = 90
        throw (1,1) double = 75
        width_HW (1,1) {mustBeNonnegative} = Inf
        width_FW (1,1) {mustBeNonnegative} = Inf
    end

    methods
        
        function self = ModelGeometry()
            % ModelGeometry Constructor to initialize the model geometry.
        end

        function [y_HW_top] = y_HW_top(self)
            % y_HW_top Returns depth of the top hanging wall, relative to mid depth.
            % Output:
            %   y_HW_top - Depth of the top hanging wall
            y_HW_top = self.get_y_top('HW');
        end

        function [y_FW_top] = y_FW_top(self)
            % y_FW_top Returns depth of the top footwall, relative to mid depth.
            % Output:
            %   y_FW_top - Depth of the top footwall
            y_FW_top = self.get_y_top('FW');
        end
 
        function [y_HW_base] = y_HW_base(self)
            % y_HW_base Returns depth of the base hanging wall, relative to mid depth.
            % Output:
            %   y_HW_base - Depth of the base hanging wall
             y_HW_base = self.get_y_base('HW');
        end

        function [y_FW_base] = y_FW_base(self)
            % y_FW_base Returns depth of the base footwall, relative to mid depth.
            % Output:
            %   y_FW_base - Depth of the base footwall
            y_FW_base = self.get_y_base('FW');
        end


        function [i_HW_top] = i_HW_top(self, y)
            % i_HW_top Returns index of the first element in the top hanging wall.
            % Input:
            %   y - Array of depths. Should be relative to depth_mid
            % Output:
            %   i_HW_top - Index of the first element in the top hanging wall
            i_HW_top = find(y <= self.y_HW_top, 1, 'first');
        end

        function [i_HW_base] = i_HW_base(self, y)
            % i_HW_base Returns index of the last element in the base hanging wall.
            % Input:
            %   y - Array of depths
            % Output:
            %   i_HW_base - Index of the last element in the base hanging wall
            i_HW_base = find(y >= self.y_HW_base, 1, 'last');
        end

        function [i_FW_top] = i_FW_top(self, y)
            % i_FW_top Returns index of the first element in the top footwall.
            % Input:
            %   y - Array of depths. Should be relative to depth_mid
            % Output:
            %   i_FW_top - Index of the first element in the top footwall
            i_FW_top = find(y <= self.y_FW_top, 1, 'first');
        end

        function [i_FW_base] = i_FW_base(self, y)
            % i_FW_base Returns index of the last element in the base footwall.
            % Input:
            %   y - Array of depths. Should be relative to depth_mid
            % Output:
            %   i_FW_base - Index of the last element in the base footwall
            i_FW_base = find(y >= self.y_FW_base, 1, 'last');
        end        
       
        function [i_FW] = i_FW(self, y)
            % i_FW Returns all indices that are wihtin the footwall interval.
            % Input:
            %   y - Array of depths. should be normalized to depth_mid
            % Output:
            %   i_FW - Indices of elements in the footwall
            i_FW = (y <= self.y_FW_top & y >= self.y_FW_base);
        end

        function [i_HW] = i_HW(self, y)
            % i_HW Returns all indices that are wihtin the hanging wall interval.
            % Input:
            %   y - Array of depths. should be normalized to depth_mid
            % Output:
            %   i_HW - Indices of elements in the hangin wall, given a
            %   certain y
            i_HW = (y <= self.y_HW_top & y >= self.y_HW_base);
        end   

        function [i_reservoir] = i_reservoir(self, y)
            % i_reservoir Returns all indices that are wihtin the hanging wall or the footwall .
            % Input:
            %   y - Array of depths. should be relative to depth_mid
            % Output:
            %   i_reservoir - Indices of elements in the hangin wall, given a
            %   certain y
            i_HW = self.i_HW(y);
            i_FW = self.i_FW(y);
            i_reservoir = i_HW | i_FW;
        end

        function [y_base] = get_y_base(self, side)
            % get_y_base Computes the base depth of the FW or HW compartment.
            % Input:
            %   side - String indicating the side ('FW' for footwall, 'HW' for hanging wall)
            % Output:
            %   y_base - Base depth of the specified compartment
            if strcmp(side, 'FW')
                y_base = -(self.thick - self.throw)/2;
            elseif strcmp(side, 'HW')
                y_base = -(self.thick + self.throw)/2;
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end
        end

        function [y_top] = get_y_top(self, side)
            % get_y_top Computes the top depth of the FW or HW compartment.
            % Input:
            %   side - String indicating the side ('FW' for footwall, 'HW' for hanging wall)
            % Output:
            %   y_top - Top depth of the specified compartment
            if strcmp(side, 'FW')
                y_top = (self.thick + self.throw)/2;
            elseif strcmp(side, 'HW')
                y_top = (self.thick - self.throw)/2; 
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end
        end

        function [L, dL] = get_along_fault_length(self, y)
            % Obtains the length along the fault, with 0 being the mid
            % depth (y=0)
            % Input:
            %   y - Depth values
            % Output
            %   L - Along-fault length
            %   dL  - Along-fault length spacing. Single value for
            %   uniformly spaced L, length(L) for varying L spacing
            L = y ./ sin(self.dip * pi/ 180);
            unique_dL = (uniquetol(diff(L), 0.001));
            if isscalar(unique_dL)
                dL = abs(unique_dL);
            else
                dL = abs(diff(L));
                dL = [dL; dL(end)];
            end
        end
        
     
    end

end