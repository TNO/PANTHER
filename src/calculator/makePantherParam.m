function p = makePantherParam(varargin)
% makePantherParam factory returning a plain struct for fast allocation
% Fields (in order): value, name, name_short, unit, uniform_with_depth,
% value_with_depth, stochastic, distribution, a, b
    props = {'value','name','name_short','unit','uniform_with_depth', 'value_with_depth', 'stochastic', 'distribution', 'a', 'b'};
    % defaults
    defaults = {[], 'undefined', '', 'undefined', true, nan, false, 'uniform', 0, 1};
    if nargin < 1
        error('Specify at least the value of the input parameter');
    end
    % build struct with defaults then override with provided args
    for i = 1:numel(props)
        p.(props{i}) = defaults{i};
    end
    for i = 1:nargin
        p.(props{i}) = varargin{i};
    end

end
