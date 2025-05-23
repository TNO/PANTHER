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
    %   top_HW_y - Returns depth of the top hanging wall, relative to mid depth
    %   top_FW_y - Returns depth of the top footwall, relative to mid depth
    %   base_HW_y - Returns depth of the base hanging wall, relative to mid depth
    %   base_FW_y - Returns depth of the base footwall, relative to mid depth
    %   top_HW_i - Returns index of the first element in the top hanging wall
    %   base_HW_i - Returns index of the last element in the base hanging wall
    %   top_FW_i - Returns index of the first element in the top footwall
    %   base_FW_i - Returns index of the last element in the base footwall
    %   get_base_y - Computes the base depth of the FW or HW compartment
    %   get_top_y - Computes the top depth of the FW or HW compartment

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

        function [top_HW_y] = top_HW_y(self)
            % top_HW_y Returns depth of the top hanging wall, relative to mid depth.
            % Output:
            %   top_HW_y - Depth of the top hanging wall
            top_HW_y = self.get_top_y('HW');
        end

        function [top_FW_y] = top_FW_y(self)
            % top_FW_y Returns depth of the top footwall, relative to mid depth.
            % Output:
            %   top_FW_y - Depth of the top footwall
            top_FW_y = self.get_top_y('FW');
        end
 
        function [base_HW_y] = base_HW_y(self)
            % base_HW_y Returns depth of the base hanging wall, relative to mid depth.
            % Output:
            %   base_HW_y - Depth of the base hanging wall
             base_HW_y = self.get_base_y('HW');
        end

        function [base_FW_y] = base_FW_y(self)
            % base_FW_y Returns depth of the base footwall, relative to mid depth.
            % Output:
            %   base_FW_y - Depth of the base footwall
            base_FW_y = self.get_base_y('FW');
        end


        function [top_HW_i] = top_HW_i(self, y)
            % top_HW_i Returns index of the first element in the top hanging wall.
            % Input:
            %   y - Array of depths
            % Output:
            %   top_HW_i - Index of the first element in the top hanging wall
            top_HW_i = find(y <= self.top_HW_y, 1, 'first');
        end

        function [base_HW_i] = base_HW_i(self, y)
            % base_HW_i Returns index of the last element in the base hanging wall.
            % Input:
            %   y - Array of depths
            % Output:
            %   base_HW_i - Index of the last element in the base hanging wall
            base_HW_i = find(y >= self.base_HW_y, 1, 'last');
        end

        function [top_FW_i] = top_FW_i(self, y)
            % top_FW_i Returns index of the first element in the top footwall.
            % Input:
            %   y - Array of depths
            % Output:
            %   top_FW_i - Index of the first element in the top footwall
            top_FW_i = find(y <= self.top_FW_y, 1, 'first');
        end

        function [base_FW_i] = base_FW_i(self, y)
            % base_FW_i Returns index of the last element in the base footwall.
            % Input:
            %   y - Array of depths
            % Output:
            %   base_FW_i - Index of the last element in the base footwall
            base_FW_i = find(y >= self.base_FW_y, 1, 'last');
        end        
       

        function [y_base] = get_base_y(self, side)
            % get_base_y Computes the base depth of the FW or HW compartment.
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

        function [y_top] = get_top_y(self, side)
            % get_top_y Computes the top depth of the FW or HW compartment.
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