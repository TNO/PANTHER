classdef PantherParameterList < FaultParameterList
    % PantherParameterList Deprecated alias for FaultParameterList.
    % See also FaultParameterList

    methods
        function self = PantherParameterList(varargin)
            warning('PantherParameterList:deprecated', ...
                'PantherParameterList is deprecated. Use FaultParameterList instead.');
            self = self@FaultParameterList(varargin{:});
        end
    end
end
