classdef TestStressChange < matlab.unittest.TestCase
    % TestStressChange Integration test for Panther
    % tests poro-elastic stress change and thermo-elastic stress change for
    % a 90 degree dip fault, without offset, against analytical expressions

    properties
    end
    
    methods (Test)
        function test_poroelastic_stress(testCase)
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
            run_instance.generate_ensemble();
            run_instance.y_extent = 0;
            y = run_instance.y;      % only evaluate at mid depth
            dx = 0;
            params = run_instance.ensemble{1};
            pressure = Pressure(params, run_instance.load_table, run_instance);
            temperature = Temperature(params, y, run_instance.load_table, 0, 'min');
            stress_change = FaultStressChange(length(y), size(pressure.dp_fault,2));
            stress_change = stress_change.calc_stress_changes(params, y, dx, pressure, temperature, 'P');
            gamma_h_infinite = (1 - 2*params.poisson)/(1 - params.poisson);   % poroelastic stress path parameter
            testCase.verifyEqual(stress_change.dsn(2), -gamma_h_infinite, "RelTol", 0.0001);
        end

         function test_thermoelastic_stress(testCase)
            % initialize run and simplify temperature steps
            run_instance = PantherInput();
            run_instance.load_table(3:end, :) = [];
            run_instance.load_table.time_steps(2) = 1;
            run_instance.load_table.T_steps(2) = -1;
            % set 0 throw, 90 degree dip
            run_instance.input_parameters.throw.value = 0;
            run_instance.input_parameters.dip.value = 90;
            run_instance.input_parameters.poisson.value = 0.2;
            run_instance.input_parameters.biot.value = 1;
            run_instance.y_extent = 0;  %  % only evaluate at mid depth
            run_instance.generate_ensemble();
            y = run_instance.y;     
            dx = 0;
            params = run_instance.ensemble{1};
            pressure = Pressure(params, run_instance.load_table, run_instance);
            temperature = Temperature(params, y, run_instance.load_table, 0, 'min');
            stress_change = FaultStressChange(length(y), size(pressure.dp_fault,2));
            stress_change = stress_change.calc_stress_changes(params, y, dx, pressure, temperature, 'T');
            gamma_h_infinite = params.young * params.therm_exp / ( 1 - params.poisson);
            testCase.verifyEqual(stress_change.dsn(2), -gamma_h_infinite, "RelTol", 0.0001);
         end

          function test_poro_thermoelastic_stress(testCase)
            % initialize run and simplify temperature steps
            run_instance = PantherInput();
            run_instance.load_table(3:end, :) = [];
            run_instance.load_table.time_steps(2) = 1;
            run_instance.load_table.T_steps(2) = -1;
            run_instance.load_table.P_steps(2) = -1;
            % set 0 throw, 90 degree dip
            run_instance.input_parameters.throw.value = 0;
            run_instance.input_parameters.dip.value = 90;
            run_instance.input_parameters.poisson.value = 0.2;
            run_instance.input_parameters.biot.value = 1;
            run_instance.y_extent = 0;  %  % only evaluate at mid depth
            run_instance.generate_ensemble();
            y = run_instance.y;      
            dx = 0;
            params = run_instance.ensemble{1};
            pressure = Pressure(params, run_instance.load_table, run_instance);
            temperature = Temperature(params, y, run_instance.load_table, 0, 'min');
            stress_change = FaultStressChange(length(y), size(pressure.dp_fault,2));
            stress_change = stress_change.calc_stress_changes(params, y, dx, pressure, temperature, 'PT');
              gamma_h_infinite_P = (1 - 2*params.poisson)/(1 - params.poisson);
            gamma_h_infinite_T = params.young * params.therm_exp / ( 1 - params.poisson);
            gamma_h_infinite = gamma_h_infinite_P + gamma_h_infinite_T;
            testCase.verifyEqual(stress_change.dsn(2), -gamma_h_infinite, "RelTol", 0.0001);
         end

         function test_thermoelastic_stress_diffusion(testCase)
            % initialize run and simplify temperature steps
            run_instance = PantherInput();
            % set 0 throw, 90 degree dip
            run_instance.input_parameters.throw.value = 0;
            run_instance.input_parameters.dip.value = 90;
            run_instance.input_parameters.poisson.value = 0.2;
            run_instance.input_parameters.biot.value = 1;
            run_instance.generate_ensemble();
            y = run_instance.y;      % only evaluate at mid depth
            dx = 0;
            params = run_instance.ensemble{1};
            pressure = Pressure(params, run_instance.load_table, run_instance);
            temperature = Temperature(params, y, run_instance.load_table, 1, 'min');
            stress_change = FaultStressChange(length(y), size(temperature.dT_fault,2));
            stress_change = stress_change.calc_stress_changes(params, y, dx, pressure, temperature, 'T');
            
         end

    end
end

