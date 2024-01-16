classdef TestPanther < matlab.unittest.TestCase
    % TestPanther Functional test for Panther

    properties
    end
    
    methods (Test)
        function test_default_single_run_P (testCase)
            result = panther();         % default inputs pressure change
            i_mid = ceil(length(result.y)/2);
            actual = result.stress{1}.sne(i_mid, end);
            expected = 30.48;
            testCase.verifyEqual(actual, expected , "RelTol", 0.1);
            actual = result.stress{1}.tau(i_mid, end);
            expected = 18.29;
            testCase.verifyEqual(actual, expected , "RelTol", 0.1);
        end

         function test_default_single_run_T (testCase)
            % test with default input, with dT and T diffusion    
            run_instance = PantherInput;
            run_instance.load_case = 'T';
            run_instance.diffusion_T = 1;
            run_instance.generate_ensemble();
            run_instance.ensemble{1}.get_gamma_T
            result = panther(run_instance);         
            i_mid = ceil(length(result.y)/2);
            actual = result.stress{1}.sne(i_mid, end);
            expected = 14.21;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            actual = result.stress{1}.tau(i_mid, end);
            expected = 7.82;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

         function test_single_with_depth_dependent_shsv(testCase)
             % test depth-variable initial stress ratio shsv
            run_instance = PantherInput;
            % set shsv varying with depth
            shsv_default = run_instance.input_parameters.shsv.value;
            run_instance.input_parameters.shsv.value_with_depth = ones(size(run_instance.y))*shsv_default;
            i_mid = ceil(length(run_instance.y)/2);
            run_instance.input_parameters.shsv.value_with_depth(i_mid - 10:i_mid + 10) = 0.8; 
            % ensure the property is set to depth-dependent (uniform = 0)
            run_instance.input_parameters.shsv.uniform_with_depth = 0;
            % run the model
            result = panther(run_instance);         
            actual = result.stress{1}.sne(i_mid, end);
            expected = 33.32;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            actual = result.stress{1}.tau(i_mid, end);
            expected = 18.62;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

    end
end

%https://github.com/marketplace/actions/run-matlab-tests