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

         function test_single_with_depth_varying_friction(testCase)
             % test depth-variable friction
            run_instance = PantherInput;
            run_instance.generate_ensemble();
            result = panther(run_instance);
            nuc_dp_uniform = result.summary.nucleation_dp;
            % make an array of f_s of size (y)
            run_instance.input_parameters.f_s.value_with_depth = ones(size(run_instance.y))*0.6;
            % set a different friction at the top of the reservoir
            run_instance.generate_ensemble();
            i_reservoir_top = run_instance.ensemble{1}.top_HW_i(run_instance.y);
            run_instance.input_parameters.f_s.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.55; 
            % ensure the property is set to depth-dependent (uniform = 0)
            run_instance.input_parameters.f_s.uniform_with_depth = 0;
            % regenerate ensemble with updated depth-varying friction
            run_instance.generate_ensemble();
            % run the model
            result = panther(run_instance);         
            actual = result.summary.nucleation_dp;
            expected = -16.5264;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            % reset to uniform friction, but with f_s of length(y)
            run_instance.input_parameters.f_s.value_with_depth = ones(size(run_instance.y))*0.6;
            run_instance.generate_ensemble();
            result = panther(run_instance);         
            actual = result.summary.nucleation_dp;
            expected = nuc_dp_uniform;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);

            % test with different f_d and d_c
            run_instance.input_parameters.f_s.value_with_depth = ones(size(run_instance.y))*0.6;
            run_instance.input_parameters.f_d.value_with_depth = ones(size(run_instance.y))*0.45;
            run_instance.input_parameters.d_c.value_with_depth = ones(size(run_instance.y))*0.005;
            % set a different friction at the top of the reservoir
            run_instance.generate_ensemble();
            i_reservoir_top = run_instance.ensemble{1}.top_HW_i(run_instance.y);
            run_instance.input_parameters.f_s.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.55;
            run_instance.input_parameters.f_d.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.43; 
            run_instance.input_parameters.d_c.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.12; 
            % ensure the property is set to depth-dependent (uniform = 0)
            run_instance.input_parameters.f_s.uniform_with_depth = 0;
            run_instance.input_parameters.f_d.uniform_with_depth = 0;
            run_instance.input_parameters.d_c.uniform_with_depth = 0;
            % regenerate ensemble with updated depth-varying friction
            run_instance.generate_ensemble();
            % run the model
            result = panther(run_instance); 
            actual = result.summary.nucleation_dp;
            expected = -18.5298;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

        function test_stochastic(testCase)
             % test depth-variable initial stress ratio shsv
            stochastic_run = PantherInput;
            % set shsv as a stochastic parameter
            stochastic_run.input_parameters.shsv.stochastic = 1;  % make shsv a stochastic parameter
            stochastic_run.input_parameters.shsv.a = 0.69;        % lower value of uniform distribution
            stochastic_run.input_parameters.shsv.b = 0.76;        % upper value of uniform distribution
            stochastic_run.stochastic = 1;                        % set the analysis to stochastic    
            stochastic_run.n_stochastic = 9;                     % number of stochastic runs
            stochastic_run.save_stress = {'first_last'};
            stochastic_run.generate_ensemble();                   % generate the ensemble (run_instance.ensemble)
            stochastic_results = panther(stochastic_run);

            actual = length(stochastic_results.stress);
            expected = 9;
            testCase.verifyEqual(actual, expected);
         end

    end
end

%https://github.com/marketplace/actions/run-matlab-tests