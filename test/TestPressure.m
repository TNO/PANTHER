classdef TestPressure < matlab.unittest.TestCase
    % test for application of fault pressure

    properties
    end
    
    methods (Test)

        function test_default_T (testCase)
            % default case
            test_case = PantherInput;
            test_case.generate_ensemble();
            pressure = PantherPressure(test_case.ensemble{1}, test_case.y,...
                test_case.load_table, test_case.load_case, test_case.diffusion_P,...
                test_case.p_fault, test_case.p_res_mode);
            % 
            i_mid.seal_seal = find(y > max(test_case.ensemble{1}.top_HW, test_case.ensemble{1}.top_FW));
            i_mid.res_base = floor((test_case.ensemble{1}.base_FW_i(y) + test_case.ensemble{1}.base_HW_i(y)) /2);

          
            % test with default input, with dT and T diffusion    
            trun = PantherInput;
            trun.load_case = 'T';
            trun.diffusion_T = 0;
            trun.generate_ensemble();
            temperature = Temperature(trun.ensemble{1}, trun.y, trun.load_table, trun.diffusion_T, 'min');
            i_mid = floor(length(trun.y)/2);
            testCase.verifyEqual(temperature.dT_fault(i_mid,end), trun.load_table.T_steps(end), "RelTol", 0.01);

            % T with linear gradient added in reservoir
            trun.input_parameters.dT_dy_multiplier.value = 0.03;
            trun.generate_ensemble();
            temperature = Temperature(trun.ensemble{1}, trun.y, trun.load_table, trun.diffusion_T, 'min');
            y_base_HW = trun.ensemble{1}.base_HW_y;
            i_base_HW = find(trun.y >= y_base_HW, 1,'last');
            add_dT = (trun.y(i_base_HW) * trun.input_parameters.dT_dy_multiplier.value);
            expected =  trun.load_table.T_steps(end) + add_dT;
            actual = temperature.dT_HW(i_base_HW, end);
            testCase.verifyEqual(actual, expected, "RelTol", 0.01);  
        end
        
%         function [i_mid] = get_mid_indices_juxtapositions(member, y)
%             
%         end
    end
end

%https://github.com/marketplace/actions/run-matlab-tests