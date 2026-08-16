classdef (HandleCompatible) FaultMesh < handle

    properties
        dy (1,1) double {mustBePositive} = 2
        y_extent (1,1) double {mustBeNonnegative} = 500
    end

    properties (Dependent)
        y (:,1) double
        faultLen (1,1) double
    end

    methods

        function [self] = FaultMesh(mesh_settings)
            if nargin == 1
                self = self.updatePropertiesFromClass(mesh_settings);
            end
        end

        function self = updatePropertiesFromClass(self, inputClass)
            % Get the list of writable stored properties of the class.
            propertyList = metaclass(self).PropertyList;
            props = {propertyList(~[propertyList.Dependent] & ~[propertyList.Constant]).Name};
            
            % Get the list of column names from the input data.
            if isstruct(inputClass)
                property_names_of_input_class = fieldnames(inputClass);
            elseif istable(inputClass)
                property_names_of_input_class = inputClass.Properties.VariableNames;
            else
                property_names_of_input_class = properties(inputClass);
            end
            
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
            ny = 1 + 2*floor(self.y_extent/self.dy);
            a = -linspace(-self.dy*fix(ny/2), self.dy*fix(ny/2), ny)';
        end

        function faultLen = get.faultLen(self)
            % faultLen returns the number of sampled y values.
            faultLen = numel(self.y);
        end

    end

end