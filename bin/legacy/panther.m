function analysis = panther(analysis)
% panther  Deprecated. Use FaultAnalyzer directly instead.
%
%   panther() duplicated what FaultAnalyzer.run() now does internally.
%   Replace any call to panther() with:
%
%       fa = FaultAnalyzer();
%       fa.run();
%
%   See also FaultAnalyzer, FaultAnalyzer.run

warning('panther:deprecated', ...
    ['panther() is deprecated and will be removed in a future release.\n' ...
     'Use FaultAnalyzer directly:\n' ...
     '    fa = FaultAnalyzer();\n' ...
     '    fa.run();']);

if nargin == 0
    analysis = FaultAnalyzer();
end
analysis = analysis.run();

end
