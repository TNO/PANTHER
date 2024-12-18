classdef ModelGeometry

    properties
        depth_mid (1,1) double {mustBeNegative} = -3000
        thick (1,1) double {mustBePositive} = 150
        dip (1,1) double = 66
        dip_azi (1,1) double = 0
        throw (1,1) double = 75
        width_HW (1,1) = Inf
        width_FW (1,1) = Inf
    end

    methods
        
        function self = ModelGeometry()
        end

        function [top_HW_y] = top_HW_y(self)
            % top_HW_y returns depth top hanging wall, relative to mid depth
            % top_HW_y = (self.thick- self.throw)/2;
            top_HW_y = get_top_y(self.thick, self.throw, 'HW');
        end

        function [top_FW_y] = top_FW_y(self)
            % top_FW_y returns depth top footwall, relative to mid depth
            % top_FW_y = (self.thick + self.throw)/2;
            top_FW_y = get_top_y(self.thick, self.throw, 'FW');
        end
 
        function [base_HW_y] = base_HW_y(self)
            % base_HW_y returns depth base hangingwall, relative to mid depth
             base_HW_y = get_base_y(self.thick, self.throw, 'HW');
        end

        function [base_FW_y] = base_FW_y(self)
            % base_FW_y returns depth base footwall, relative to mid depth
            % base_FW_y = -(self.thick - self.throw)/2;
            base_FW_y = get_base_y(self.thick, self.throw, 'FW');
        end


        function [top_HW_i] = top_HW_i(self, y)
            % top_HW_y returns index of first element in top hanging wall 
            top_HW_i = find(y <= self.top_HW_y, 1, 'first');
        end

        function [base_HW_i] = base_HW_i(self, y)
            % top_FW_i index base footwall 
            base_HW_i = find(y >= self.base_HW_y, 1, 'last');
        end

        function [top_FW_i] = top_FW_i(self, y)
            % top_HW_y returns index of first element in top hanging wall 
            top_FW_i = find(y <= self.top_FW_y, 1, 'first');
        end

        function [base_FW_i] = base_FW_i(self, y)
            % top_FW_i index base footwall 
            base_FW_i = find(y >= self.base_FW_y, 1, 'last');
        end        
       

        function [y_base] = get_base_y(self, side)
            % function to compute the base depth of the FW or HW compartment, with
            % resepct to the center depth (depth_ymid)
            % thick = [double] reservoir thickness
            % throw = [double] reservoir throw
            % side = [string] side, FW: footwall and HW: hanging wall
            if strcmp(side, 'FW')
                y_base = -(self.thick - self.throw)/2;
            elseif strcmp(side, 'HW')
                y_base = -(self.thick + self.throw)/2;
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end
        end

        function [y_top] = get_top_y(self, side)
            % function to compute the base top of the FW or HW compartment, with
            % resepct to the center depth (depth_ymid)
            % thick = [double] reservoir thickness
            % throw = [double] reservoir throw
            % side = [string] side, FW: footwall and HW: hanging wall
            if strcmp(side, 'FW')
                y_top = (self.thick + self.throw)/2;
            elseif strcmp(side, 'HW')
                y_top = (self.thick - self.throw)/2; 
            else
                error('Incorrect side indicator entered, should be FW or HW');
            end
        end
        
     
    end

end