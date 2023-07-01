classdef TestStressChange < matlab.unittest.TestCase
    % TestPanther Integration test for Panther

    properties
    end
    
    methods (Test)
        function test_poroelastic_stress(testCase)
            % initialize run and simplify pressure steps
            run_instance = PantherInput();
            run_instance.diffusion_P = 0;
            run_instance.load_table(3:end, :) = [];
            run_instance.load_table.time_steps(2) = 1;
            run_instance.load_table.P_steps(2) = -1;
            % set 0 throw, 90 degree dip
            run_instance.input_parameters.throw.value = 0;
            run_instance.input_parameters.dip.value = 90;
            run_instance.input_parameters.poisson.value = 0.2;
            run_instance.input_parameters.biot.value = 1;
            run_instance.generate_ensemble();
            y = run_instance.y;
            y = 0;
            dx = 0;
            params = run_instance.ensemble{1};
            pressure = PantherPressure(params, y, run_instance.load_table, 'P', 0, 'min');
            stress_change = FaultStressChange(length(y), size(pressure.dp_fault,2));
            stress_change = stress_change.calc_stress_changes(params, y, dx, pressure, 'P');
            gamma_h_infinite = (1 - 2*0.2)/(1 - 0.2);   % poroelastic stress path parameter
            testCase.verifyEqual(stress_change.dsn(2), gamma_h_infinite, "RelTol", 0.0001);
        end
    end
end

