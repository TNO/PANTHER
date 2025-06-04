classdef TestPressure < matlab.unittest.TestCase
    % test for application of fault pressure

    properties
    end
    
    methods (Test)

        function test_default_P (testCase)
            % default case
            tc = PantherAnalysis;
            tc.load_table = tc.load_table(1:2,:);
            tc.load_table.time_steps(2) = 1;
            tc.load_table.P_steps(2) = -1;
            tc.generate_ensemble();
            
            % 
            i_mid.seal_seal = floor(find(tc.y > max(tc.ensemble_members{1}.top_HW_y, tc.ensemble_members{1}.top_FW_y),1,'last')/2);
            i_mid.res_res = floor(length(tc.y)/2) ;
            i_mid.res_base = floor((tc.ensemble_members{1}.base_FW_i(tc.y) + tc.ensemble_members{1}.base_HW_i(tc.y)) /2);
            i_mid.res_seal = floor((tc.ensemble_members{1}.top_FW_i(tc.y) + tc.ensemble_members{1}.top_HW_i(tc.y)) /2);
            i_mid.base_base = floor(find(tc.y < min(tc.ensemble_members{1}.base_HW_y, tc.ensemble_members{1}.base_FW_y),1,'first'));
            
            % case t < h, P_res_mode = 'same', P_fault_mode = 'min', diffusion=0
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);
            testCase.verifyEqual(p.dP(i_mid.res_res, end), -1, "RelTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_base, end), -1, "RelTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_seal, end), -1, "RelTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.base_base, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.seal_seal, end), 0, "AbsTol", 1e-10);
            
            % set p fault to max(p_Hw, p_FW)
            tc.P_fault_mode = 'max';
            tc.generate_ensemble;
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);
          
            testCase.verifyEqual(p.dP(i_mid.res_res, end), -1, "RelTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_base, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_seal, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.base_base, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.seal_seal, end), 0, "AbsTol", 1e-10);
            
            % single side scenario
            tc.P_fault_mode = 'min';
            tc.input_parameters.width_FW.value = 0;
            tc.input_parameters.width_HW.value = inf;
            tc.generate_ensemble;
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);
            
            testCase.verifyEqual(p.dP(i_mid.res_res, end), -1, "RelTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_base, end), -1, "RelTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_seal, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.base_base, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.seal_seal, end), 0, "AbsTol", 1e-10);

            % single side scenario
            tc.P_fault_mode = 'min';
            tc.input_parameters.width_FW.value = inf;
            tc.input_parameters.width_HW.value = 0;
            tc.generate_ensemble;
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);
            
            testCase.verifyEqual(p.dP(i_mid.res_res, end), -1, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_base, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.res_seal, end), -1, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.base_base, end), 0, "AbsTol", 1e-10);
            testCase.verifyEqual(p.dP(i_mid.seal_seal, end), 0, "AbsTol", 1e-10);

        end

        function test_pressure_diffusion(testCase)
            % default case
            tc = PantherAnalysis;
            tc.load_table = tc.load_table(1:2,:);
            tc.load_table.time_steps(2) = 10;
            tc.load_table.P_steps(2) = -10;
            tc.diffusion_P = 1;
            tc.input_parameters.P_over.value = 2;
            tc.P_fault_mode = 'min';
            tc.generate_ensemble();
            
            % case t < h, P_res_mode = 'same', p_fault = 'min', diffusion =
            % 1, p_over = 2 MPa
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);

            y_eval = -150;
            p_at_eval = interp1(tc.y, p.P(:,end), y_eval);

            testCase.verifyEqual(p_at_eval, 31.0805, "RelTol", 0.05);
        end

        function test_initial_pressures(testCase)
            % default case
            tc = PantherAnalysis;
            tc.load_table = tc.load_table(1:2,:);
            tc.load_table.time_steps(2) = 1;
            tc.load_table.P_steps(2) = -1;
            tc.diffusion_P = 1;
            tc.input_parameters.P_grad_res.value = 0.2;
            tc.P0_fault_mode = 'max';
            tc.P_fault_mode = 'min';
            tc.generate_ensemble();
            
            % case t < h, P_res_mode = 'same', dP_fault_mode = 'min', diffusion =
            % 1, P_fault_mode = 'max'
            % initial pressure should be equal to pressure computed at t =
            % 0 (sanity checks, pressures at t=0 should be equal, dp should be 0)
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);
            P_HW = p.get_P_HW();
            P0_HW = p.get_P0_HW();
            P_HW_diff = P_HW(:,1) - P0_HW;
            testCase.verifyEqual(max(abs(P_HW_diff)), 0, "RelTol", 1e-10);
            dp_HW = p.get_dP_HW();
            testCase.verifyEqual(max(abs(dp_HW(:,1))), 0, "RelTol", 1e-10);
            dp_FW = p.get_dP_FW();
            testCase.verifyEqual(max(abs(dp_FW(:,1))), 0, "RelTol", 1e-10);
            testCase.verifyEqual(max(abs(p.dP(:,1))), 0, "RelTol", 1e-10);
        end
        
        function test_pressure_setting_for_different_p_grad(testCase)
             % default case
            tc = PantherAnalysis;
            tc.load_table = tc.load_table(1:2,:);
            tc.load_table.time_steps(2) = 1;
            tc.load_table.P_steps(2) = -1;
            tc.diffusion_P = 0;
            tc.input_parameters.P_grad_res.value = 0.2;
            tc.input_parameters.P_over.value = 0.2;
            tc.P_fault_mode = 'max';
            tc.P0_fault_mode = 'min';
            tc.P_res_mode = 'same';
            tc.generate_ensemble();
            p = Pressure(tc.ensemble_members{1}, tc.load_table, tc);
            % seal reservoir juxtaposition
            i_seal_res = floor((p.top_FW_i(p.y) + p.top_HW_i(p.y))/2);
            % expected pressure equal to p gradient without p_grad_res and
            % p_over because P_fault_mode = 'min' 
            expected = -(1/1000)*(p.y(i_seal_res) + tc.input_parameters.depth_mid.value)...,
                *tc.input_parameters.P_grad.value + tc.input_parameters.P_offset.value;
            observed = p.P(i_seal_res, 1) ;
            testCase.verifyEqual(observed, expected, "RelTol", 1e-10);
        end

    end
end

%https://github.com/marketplace/actions/run-matlab-tests