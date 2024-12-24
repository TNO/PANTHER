classdef TestPressure < matlab.unittest.TestCase
    % test for application of fault pressure

    properties
    end
    
    methods (Test)

        function test_default_P (testCase)
            % default case
            tc = PantherInput;
            tc.load_table = tc.load_table(1:2,:);
            tc.load_table.time_steps(2) = 1;
            tc.load_table.P_steps(2) = -1;
            tc.generate_ensemble();
            
            % 
            i_mid.seal_seal = floor(find(tc.y > max(tc.ensemble{1}.top_HW_y, tc.ensemble{1}.top_FW_y),1,'last')/2);
            i_mid.res_res = floor(length(tc.y)/2) ;
            i_mid.res_base = floor((tc.ensemble{1}.base_FW_i(tc.y) + tc.ensemble{1}.base_HW_i(tc.y)) /2);
            i_mid.res_seal = floor((tc.ensemble{1}.top_FW_i(tc.y) + tc.ensemble{1}.top_HW_i(tc.y)) /2);
            i_mid.base_base = floor(find(tc.y < min(tc.ensemble{1}.base_HW_y, tc.ensemble{1}.base_FW_y),1,'first'));
            
            % case t < h, p_res_mode = 'same', p_fault = 'min', diffusion=0
            p = Pressure(tc.ensemble{1}, tc.load_table, tc);
            testCase.verifyEqual(p.dp_fault(i_mid.res_res, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.res_base, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.res_seal, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.base_base, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.seal_seal, end), 0);
            
            % set dp fault to max(dp_Hw, dp_FW)
            tc.dp_fault_mode = 'max';
            tc.generate_ensemble;
            p = Pressure(tc.ensemble{1}, tc.load_table, tc);
          
            testCase.verifyEqual(p.dp_fault(i_mid.res_res, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.res_base, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.res_seal, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.base_base, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.seal_seal, end), 0);
            
            % single side scenario
            tc.dp_fault_mode = 'min';
            tc.input_parameters.width_FW.value = 0;
            tc.input_parameters.width_HW.value = inf;
            tc.generate_ensemble;
            p = Pressure(tc.ensemble{1}, tc.load_table, tc);
            
            testCase.verifyEqual(p.dp_fault(i_mid.res_res, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.res_base, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.res_seal, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.base_base, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.seal_seal, end), 0);

            % single side scenario
            tc.p_fault_mode = 'min';
            tc.input_parameters.width_FW.value = inf;
            tc.input_parameters.width_HW.value = 0;
            tc.generate_ensemble;
            p = Pressure(tc.ensemble{1}, tc.load_table, tc);
            
            testCase.verifyEqual(p.dp_fault(i_mid.res_res, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.res_base, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.res_seal, end), -1);
            testCase.verifyEqual(p.dp_fault(i_mid.base_base, end), 0);
            testCase.verifyEqual(p.dp_fault(i_mid.seal_seal, end), 0);

             % dp_HW = p.get_HW_pressure_change();
             % dp_FW = p.get_FW_pressure_change();
             % plot(dp_HW,tc.y, dp_FW, tc.y)

        end

        function test_pressure_diffusion(testCase)
            % default case
            tc = PantherInput;
            tc.load_table = tc.load_table(1:2,:);
            tc.load_table.time_steps(2) = 10;
            tc.load_table.P_steps(2) = -10;
            tc.diffusion_P = 1;
            tc.input_parameters.p_over.value = 2;
            tc.dp_fault_mode = 'max_abs';
            tc.generate_ensemble();
            
            % case t < h, p_res_mode = 'same', p_fault = 'min', diffusion =
            % 1, p_over = 2 MPa
            p = Pressure(tc.ensemble{1}, tc.load_table, tc);

            y_eval = -150;
            p_at_eval = interp1(tc.y, p.p(:,end), y_eval);

            testCase.verifyEqual(p_at_eval, 31.0805, "RelTol", 0.05);
        end
        

    end
end

%https://github.com/marketplace/actions/run-matlab-tests