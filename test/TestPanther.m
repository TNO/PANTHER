classdef TestPanther < matlab.unittest.TestCase
    % TestPanther Functional test for Panther

    properties
    end
    
    methods (Test)
        function test_default_single_run_P (testCase)
            result = PantherAnalysis();
            result = result.run();
            i_mid = ceil(length(result.y)/2);
            actual = result.faultResults.sne(i_mid, end);
            expected = 30.48;
            testCase.verifyEqual(actual, expected , "RelTol", 0.1);
            actual = result.faultResults.tau(i_mid, end);
            expected = 18.29;
            testCase.verifyEqual(actual, expected , "RelTol", 0.1);
        end

         function test_default_single_run_T (testCase)
            % test with default input, with dT and T diffusion    
            run_instance = PantherAnalysis;
            run_instance.load_case = 'T';
            run_instance.diffusion_T = 0;
            run_instance.generate_ensemble();
            run_instance.ensemble_members{1}.get_gamma_T;
                run_instance = run_instance.run();
            i_mid = ceil(length(run_instance.y)/2);
                actual = run_instance.faultResults.sne(i_mid, end);
            expected = 14.29;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
                actual = run_instance.faultResults.tau(i_mid, end);
            expected = 8.055;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

         function test_single_with_depth_dependent_shsv(testCase)
             % test depth-variable initial stress ratio shsv
            run_instance = PantherAnalysis;
            % set shsv varying with depth
            shsv_default = run_instance.getInputParameter('shsv');
            shsv_with_depth = ones(size(run_instance.y))*shsv_default;
            i_mid = ceil(length(run_instance.y)/2);
            shsv_with_depth(i_mid - 10:i_mid + 10) = 0.8;
            run_instance.setDepthDependentInputParameter('shsv', shsv_with_depth);
            % run the model
                run_instance = run_instance.run();
                actual = run_instance.faultResults.sne(i_mid, end);
            expected = 33.32;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
                actual = run_instance.faultResults.tau(i_mid, end);
            expected = 18.62;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

         function test_single_with_depth_varying_friction(testCase)
             % test depth-variable friction
            run_instance = PantherAnalysis;
            run_instance.run();
            i_mid = ceil(length(run_instance.y)/2);
            nuc_step = run_instance.faultSummary.nucleation_load_step;
            nuc_dp_uniform = run_instance.faultResults.dP(i_mid, nuc_step);
            % make an array of f_s of size (y)
            f_s_with_depth = ones(size(run_instance.y))*0.6;
            % set a different friction at the top of the reservoir
            run_instance.generate_ensemble();
            i_reservoir_top = run_instance.ensemble_members{1}.i_HW_top(run_instance.y);
            f_s_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.55;
            run_instance.setDepthDependentInputParameter('f_s', f_s_with_depth);
            % run the model
            run_instance.generate_ensemble();
            run_instance = run_instance.run();
            nuc_step = run_instance.faultSummary.nucleation_load_step;
            actual = run_instance.faultResults.dP(i_mid, nuc_step);
            expected = -17.53;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            % reset to uniform friction, but with f_s of length(y)
            run_instance.setDepthDependentInputParameter('f_s', ones(size(run_instance.y))*0.6);
            run_instance = run_instance.run();
            nuc_step = run_instance.faultSummary.nucleation_load_step;
            actual = run_instance.faultResults.dP(i_mid, nuc_step);
            expected = nuc_dp_uniform;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
            % test with different f_d and d_c
            f_s_with_depth = ones(size(run_instance.y))*0.6;
            f_d_with_depth = ones(size(run_instance.y))*0.45;
            d_c_with_depth = ones(size(run_instance.y))*0.005;
            % set a different friction at the top of the reservoir
            i_reservoir_top = run_instance.ensemble_members{1}.i_HW_top(run_instance.y);
            f_s_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.55;
            f_d_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.43; 
            d_c_with_depth(i_reservoir_top - 15: i_reservoir_top + 15) = 0.12;
            run_instance.setDepthDependentInputParameter('f_s', f_s_with_depth);
            run_instance.setDepthDependentInputParameter('f_d', f_d_with_depth);
            run_instance.setDepthDependentInputParameter('d_c', d_c_with_depth);
            % run the model
            run_instance = run_instance.run();
            nuc_step = run_instance.faultSummary.nucleation_load_step;
            actual = run_instance.faultResults.dP(i_mid, nuc_step);
            % expected = -19.77;
            % expected = -20.50;
            expected = -21.04;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);
         end

        function test_sH_dir(testCase)
            % test to check whether sH_dir is handled correctly
            % initialize run and simplify pressure steps
            run_instance = PantherAnalysis();
            run_instance.load_table(3:end, :) = [];
            run_instance.load_table.time_steps(2) = 1;
            run_instance.load_table.P_steps(2) = -1;
            % set 0 throw, 90 degree dip
            run_instance.setInputParameter('throw', 0);
            run_instance.setInputParameter('dip', 90);
            run_instance.setInputParameter('poisson', 0.2);
            run_instance.setInputParameter('biot', 1);
            run_instance.setInputParameter('sv_grad', 22);
            run_instance.setInputParameter('shsv', 0.75);
            run_instance.setInputParameter('sHsh', 1.1);
            run_instance.y_extent = 0;

            % test with strike parallel to sH_dir
            run_instance.setInputParameter('sH_dir', 0);
            run_instance.setInputParameter('dip_azi', 90);   % strike parallel to sH_dir
            result = run_instance.run();
            actual = result.faultResults.sne(1);
            expected = -run_instance.ensemble_members{1}.depth_mid/1000 * ((run_instance.ensemble_members{1}.sv_grad ...
                * run_instance.ensemble_members{1}.shsv) - run_instance.ensemble_members{1}.P_grad) ; 
            testCase.verifyEqual(actual, expected, "RelTol", 1e-10);
            
            % test with sH_dir perpendicular to strike (parallel to
            % dip_azi)
            run_instance.setInputParameter('sH_dir', 90);
            result = run_instance.run();
            actual = result.faultResults.sne(1);
            expected = -run_instance.ensemble_members{1}.depth_mid/1000 * ((run_instance.ensemble_members{1}.sv_grad ...
                * run_instance.ensemble_members{1}.shsv * run_instance.ensemble_members{1}.sHsh) - run_instance.ensemble_members{1}.P_grad) ; 
            testCase.verifyEqual(actual, expected, "RelTol", 1e-10); 

            % test with sH_dir perpendicular to strike (parallel to
            % dip_azi), with sH_dir given as negative
            run_instance.setInputParameter('sH_dir', -90);
            result = run_instance.run();
            actual = result.faultResults.sne(1);
            expected = -run_instance.ensemble_members{1}.depth_mid/1000 * ((run_instance.ensemble_members{1}.sv_grad ...
                * run_instance.ensemble_members{1}.shsv * run_instance.ensemble_members{1}.sHsh) - run_instance.ensemble_members{1}.P_grad) ; 
            testCase.verifyEqual(actual, expected, "RelTol", 1e-10); 
        end

    end
end