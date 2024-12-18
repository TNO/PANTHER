classdef FaultMesh

    properties
        dy (1,1) double {mustBePositive} = 2
        y_extent (1,1) double {mustBePositive} = 500
    end

    properties (Dependent)
        y
    end

    methods

        function [self] = FaultMesh(mesh_settings)
            if nargin == 1
                self = self.updatePropertiesFromClass(mesh_settings);
            end
        end

        function self = updatePropertiesFromClass(self, inputClass)
            % Get the list of properties of the class
            props = properties(self);
            
            % Get the list of column names from the input table
            property_names_of_input_class = properties(inputClass);
            
            % Loop through each property and update if there's a matching column in the table
            for i = 1:length(props)
                property_name = props{i};
                if ismember(property_name, property_names_of_input_class)
                    self.(property_name) = inputClass.(property_name);
                end
            end
        end

        function a = get.y(self)
            % y initializes depth y relative to depth_mid
            a = get_fault_y(self.dy, self.y_extent);
        end

        function y = get_fault_y(~, dy, y_extent)
            ny = 1 + 2*floor(y_extent/dy);
            dy_2 = 2*y_extent/ny;           % Rescale dy
            y = zeros(ny, 1);
            
            for i = 1 : ny
                y(i) = y_extent - (i-1)*dy_2;
            end
        end

    end

end