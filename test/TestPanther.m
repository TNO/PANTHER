classdef TestPanther < matlab.unittest.TestCase
    % TestPanther Integration test for Panther

    properties
    end
    
    methods (Test)
        function test_default_single_run (testCase)
            result = panther();         % default inputs
            i_mid = ceil(length(result.y)/2);
            actual = result.stress{1}.sne(i_mid, end);
            expected = 26.668;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            actual = result.stress{1}.tau(i_mid, end);
            expected = 16.0;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
        end
    end
end

%https://github.com/marketplace/actions/run-matlab-tests