classdef FaultSurfaceCalculator
% FaultSurfaceCalculator  Deprecated. Use MultiFaultAnalyzer instead.
%
%   FaultSurfaceCalculator has been superseded by MultiFaultAnalyzer, which
%   provides the same functionality with a cleaner, more consistent API.
%
%   Replace:
%       fsc = FaultSurfaceCalculator(n);
%       fsc.run();
%
%   With:
%       mfa = MultiFaultAnalyzer();
%       mfa.initialize(n);
%       mfa.run();
%
%   See also MultiFaultAnalyzer

    methods
        function self = FaultSurfaceCalculator(varargin) %#ok<VANUS>
            warning('FaultSurfaceCalculator:deprecated', ...
                ['FaultSurfaceCalculator is deprecated and will be removed in a future release.\n' ...
                 'Use MultiFaultAnalyzer instead:\n' ...
                 '    mfa = MultiFaultAnalyzer();\n' ...
                 '    mfa.initialize(n);']);
        end
    end

end
