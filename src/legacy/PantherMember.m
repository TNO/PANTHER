classdef PantherMember < FaultRealization
    % PantherMember  Deprecated alias for FaultRealization.
    % See also FaultRealization
    methods
        function self = PantherMember(varargin)
            warning('PantherMember:deprecated', ...
                'PantherMember is deprecated. Use FaultRealization instead.');
            self = self@FaultRealization(varargin{:});
        end
    end
end
