classdef PantherParam
    % PantherParam create Panther input parameter object

    properties
        value double                % parameter default value
        name string {mustBeText}    % parameter name e.g. 'Young's modulus'
        unit string {mustBeText}    % parameter unit e.g. 'Pa'       
        uniform_with_depth logical = 1;   % uniform with depth (NB non-uniform not used yet)
        value_with_depth double     % 
        stochastic logical = 0;     % stochastic analysis drawing from distribution
        distribution {mustBeMember(distribution,{'uniform','logunifom'})} = 'uniform'; % distribution from which to sample values for the stochastic analysis
        a double = 0;               % first parameter used to create distribution. e.g. for uniform, low
        b double = 1;               % second parameter used to create distribution. e.g. for uniform, high
    end

    methods 
        function self = PantherParam(varargin)
            props = properties(self);
            % nargin TODO check at least 3 inputs are given
            for i = 1 : nargin  
                self.(props{i}) = varargin{i};
            end
        end
        function self = SetStochastic(self, dist, a, b)
            self.distribution = dist;
            self.a = a;
            self.b = b;
        end
        % TODO add simple depth-dependency for non-uniform
    end

end