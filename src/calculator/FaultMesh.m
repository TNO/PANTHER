classdef (HandleCompatible) FaultMesh < handle

    properties
        dy (1,1) double {mustBePositive} = 2
        y_extent (1,1) double {mustBeNonnegative} = 500
    end

    properties (Dependent, Hidden)
        y (:,1) double
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
            a = self.get_fault_y(self.dy, self.y_extent);
        end

        function y = get_fault_y(~, dy, y_extent)
            % Computes the depth array y
            % y will run to the nearest value specified by y_extent
            % (dy will not be rescaled to precisely match y_extent)
            % Input
            %   dy:         depth spacing
            %   y_extent:   depth extent +- with respect to zero
            % Output
            %   y:          depth array. runs from +ve (upwards) to -ve
            %   (downwards) 
            ny = 1 + 2*floor(y_extent/dy);
            y = -linspace(-dy*fix(ny/2), dy*fix(ny/2), ny)';
        end

    end

end