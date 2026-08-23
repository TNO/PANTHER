classdef (HandleCompatible) PantherAnalysis < FaultAnalyzer
    % PantherAnalysis  Deprecated alias for FaultAnalyzer.
    %
    % This class exists solely for backward compatibility. All new code
    % should use FaultAnalyzer directly.
    %
    % See also FaultAnalyzer

    methods
        function self = PantherAnalysis(varargin)
            % Issue a one-time warning per session and forward to FaultAnalyzer.
            warning('PantherAnalysis:deprecated', ...
                ['PantherAnalysis is deprecated and will be removed in a future release. ' ...
                 'Use FaultAnalyzer instead.']);
            self = self@FaultAnalyzer(varargin{:});
        end
    end

end
