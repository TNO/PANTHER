classdef TestTemperature < matlab.unittest.TestCase
    % integration tests for temperature

    properties
    end
    
    methods (Test)

        function test_default_T (testCase)
            % test with default input, with dT and T diffusion    
            trun = PantherAnalysis;
            trun.load_case = 'T';
            trun.diffusion_T = 0;
            trun.generate_ensemble();
            temperature = Temperature(trun, 'min');
            i_mid = floor(length(trun.y)/2);
            testCase.verifyEqual(temperature.dT(i_mid,end), trun.load_table.T_steps(end), "RelTol", 0.01);

            % T with linear gradient added in reservoir
            trun.setInputParameter('dT_dy_multiplier', 0.03);
            trun.generate_ensemble();
            temperature = Temperature(trun, 'min');
            y_base_HW = trun.ensemble_members{1}.y_HW_base;
            i_base_HW = find(trun.y >= y_base_HW, 1,'last');
            add_dT = (trun.y(i_base_HW) * trun.getInputParameter('dT_dy_multiplier'));
            expected =  trun.load_table.T_steps(end) + add_dT;
            actual = temperature.dT_HW(i_base_HW, end);
            testCase.verifyEqual(actual, expected, "RelTol", 0.01);  
        end
    end
end

%https://github.com/marketplace/actions/run-matlab-tests