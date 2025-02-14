classdef PantherParam
    % PantherParam create Panther input parameter object

    properties
        value double                        % parameter default value
        name string {mustBeText} = 'undefined'   % parameter name e.g. 'Young's modulus'
        unit string {mustBeText} = 'undefined'   % parameter unit e.g. 'Pa'       
        uniform_with_depth logical = 1;     % uniform with depth
        value_with_depth double             % 
        stochastic logical = 0;             % stochastic analysis drawing from distribution
        distribution {mustBeMember(distribution,{'uniform','logunifom'})} = 'uniform'; % distribution from which to sample values for the stochastic analysis
        a double = 0;                       % first parameter used to create distribution. e.g. for uniform, low
        b double = 1;                       % second parameter used to create distribution. e.g. for uniform, high
    end

    methods 
        function self = PantherParam(varargin)
            % constructor for Parameter object
            props = properties(self);
            disp(length(varargin))
            if length(varargin) < 1
                error('Specify at least the value of the input parameter');
            end
            for i = 1 : nargin  
                self.(props{i}) = varargin{i};
            end
        end

        function self = SetStochastic(self, dist, a, b)
            self.distribution = dist;
            self.a = a;
            self.b = b;
        end

    end

end