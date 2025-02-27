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
            run_instance.diffusion_T = 0;
            run_instance.generate_ensemble();
            run_instance.ensemble_members{1}.get_gamma_T;
            run_instance = panther(run_instance);         
            i_mid = ceil(length(run_instance.y)/2);
            actual = run_instance.stress{1}.sne(i_mid, end);
            expected = 14.29;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            actual = run_instance.stress{1}.tau(i_mid, end);
            expected = 8.055;
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
            run_instance = panther(run_instance);         
            actual = run_instance.stress{1}.sne(i_mid, end);
            expected = 33.32;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            actual = run_instance.stress{1}.tau(i_mid, end);
            expected = 18.62;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

         function test_single_with_depth_varying_friction(testCase)
             % test depth-variable friction
            run_instance = PantherInput;
            run_instance = panther(run_instance);
            nuc_dp_uniform = run_instance.summary.nucleation_dP;
            % make an array of f_s of size (y)
            run_instance.input_parameters.f_s.value_with_depth = ones(size(run_instance.y))*0.6;
            % set a different friction at the top of the reservoir
            run_instance.generate_ensemble();
            i_reservoir_top = run_instance.ensemble_members{1}.top_HW_i(run_instance.y);
            run_instance.input_parameters.f_s.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.55;
            % ensure the property is set to depth-dependent (uniform = 0)
            run_instance.input_parameters.f_s.uniform_with_depth = 0;
            % run the model
            run_instance.generate_ensemble()
            run_instance = panther(run_instance);         
            actual = run_instance.summary.nucleation_dP;
            expected = -17.53;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            % reset to uniform friction, but with f_s of length(y)
            run_instance.input_parameters.f_s.value_with_depth = ones(size(run_instance.y))*0.6;
            run_instance = panther(run_instance);         
            actual = run_instance.summary.nucleation_dP;
            expected = nuc_dp_uniform;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            % test with different f_d and d_c
            run_instance.input_parameters.f_s.value_with_depth = ones(size(run_instance.y))*0.6;
            run_instance.input_parameters.f_d.value_with_depth = ones(size(run_instance.y))*0.45;
            run_instance.input_parameters.d_c.value_with_depth = ones(size(run_instance.y))*0.005;
            % set a different friction at the top of the reservoir
            i_reservoir_top = run_instance.ensemble_members{1}.top_HW_i(run_instance.y);
            run_instance.input_parameters.f_s.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.55;
            run_instance.input_parameters.f_d.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.43; 
            run_instance.input_parameters.d_c.value_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.12; 
            % ensure the property is set to depth-dependent (uniform = 0)
            run_instance.input_parameters.f_s.uniform_with_depth = 0;
            run_instance.input_parameters.f_d.uniform_with_depth = 0;
            run_instance.input_parameters.d_c.uniform_with_depth = 0;
            % run the model
            run_instance = panther(run_instance); 
            actual = run_instance.summary.nucleation_dP;
            %expected = -19.77;
            expected = -20.80;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

        function test_stochastic(testCase)
             % test depth-variable initial stress ratio shsv
            stochastic_run = PantherInput();
            % set shsv as a stochastic parameter
            stochastic_run.input_parameters.shsv.stochastic = 1;  % make shsv a stochastic parameter
            stochastic_run.input_parameters.shsv.a = 0.69;        % lower value of uniform distribution
            stochastic_run.input_parameters.shsv.b = 0.76;        % upper value of uniform distribution
            stochastic_run.stochastic = 1;                        % set the analysis to stochastic    
            stochastic_run.n_stochastic = 9;                     % number of stochastic runs
            stochastic_run.save_stress = {'first_last'};
            stochastic_run = panther(stochastic_run);

            actual = length(stochastic_run.stress);
            expected = 9;
            testCase.verifyEqual(actual, expected);
         end

        function test_sH_dir(testCase)
            % test to check whether sH_dir is handled correctly
            % initialize run and simplify pressure steps
            run_instance = PantherInput();
            run_instance.load_table(3:end, :) = [];
            run_instance.load_table.time_steps(2) = 1;
            run_instance.load_table.P_steps(2) = -1;
            % set 0 throw, 90 degree dip
            run_instance.input_parameters.throw.value = 0;
            run_instance.input_parameters.dip.value = 90;
            run_instance.input_parameters.poisson.value = 0.2;
            run_instance.input_parameters.biot.value = 1;
            run_instance.input_parameters.sv_grad.value = 22;
            run_instance.input_parameters.shsv.value = 0.75;
            run_instance.input_parameters.sHsh.value = 1.1;
            run_instance.y_extent = 0;

            % test with strike parallel to sH_dir
            run_instance.input_parameters.sH_dir.value = 0;
            run_instance.input_parameters.dip_azi.value = 90;   % strike parallel to sH_dir
            result = panther(run_instance);
            actual = result.stress{1}.sne(1); 
            expected = -run_instance.ensemble_members{1}.depth_mid/1000 * ((run_instance.ensemble_members{1}.sv_grad ...
                * run_instance.ensemble_members{1}.shsv) - run_instance.ensemble_members{1}.P_grad) ; 
            testCase.verifyEqual(actual, expected, "RelTol", 1e-10);
            
            % test with sH_dir perpendicular to strike (parallel to
            % dip_azi)
            run_instance.input_parameters.sH_dir.value = 90;
            result = panther(run_instance);
            actual = result.stress{1}.sne(1);
            expected = -run_instance.ensemble_members{1}.depth_mid/1000 * ((run_instance.ensemble_members{1}.sv_grad ...
                * run_instance.ensemble_members{1}.shsv * run_instance.ensemble_members{1}.sHsh) - run_instance.ensemble_members{1}.P_grad) ; 
            testCase.verifyEqual(actual, expected, "RelTol", 1e-10); 

            % test with sH_dir perpendicular to strike (parallel to
            % dip_azi), with sH_dir given as negative
            run_instance.input_parameters.sH_dir.value = -90;
            result = panther(run_instance);
            actual = result.stress{1}.sne(1);
            expected = -run_instance.ensemble_members{1}.depth_mid/1000 * ((run_instance.ensemble_members{1}.sv_grad ...
                * run_instance.ensemble_members{1}.shsv * run_instance.ensemble_members{1}.sHsh) - run_instance.ensemble_members{1}.P_grad) ; 
            testCase.verifyEqual(actual, expected, "RelTol", 1e-10); 
        end

    end
end