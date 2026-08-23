function mfa = MultiFaultCalculator(n_pillars)
% MultiFaultCalculator  Deprecated. Use MultiFaultAnalyzer instead.
%
%   Replace:
%       mfc = MultiFaultCalculator(n);
%       mfc.run();
%
%   With:
%       mfa = MultiFaultAnalyzer();
%       mfa.initialize(n);
%       mfa.run();
%
%   See also MultiFaultAnalyzer

warning('MultiFaultCalculator:deprecated', ...
    ['MultiFaultCalculator is deprecated and will be removed in a future release.\n' ...
     'Use MultiFaultAnalyzer instead:\n' ...
     '    mfa = MultiFaultAnalyzer();\n' ...
     '    mfa.initialize(n);']);

mfa = MultiFaultAnalyzer();
if nargin > 0
    mfa.initialize(n_pillars);
end

end
