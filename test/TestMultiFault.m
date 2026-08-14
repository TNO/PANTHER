classdef TestMultiFault < matlab.unittest.TestCase
    % TestMultiFault Functional test for MultiFault controller

    properties
    end
    
    methods (Test)
        function test_input_parameter_getters_and_setters(testCase)
            n_faults = 100;
            mf = MultiFault();
            mf = mf.initialize(n_faults);
            default_dips = mf.get_input_parameter('dip');
            testCase.verifyEqual(length(default_dips), n_faults);
            % ensure the following do not throw errors
            mf.faults(2).set_input_parameter('dip', 80);
            testCase.verifyEqual(max(mf.get_input_parameter('dip')), 80);
            mf.set_input_parameter('dip',default_dips);
            testCase.verifyEqual(max(mf.get_input_parameter('dip')), 70);
        end

        function test_depth_dependent_input_parameter_getters_and_setters(testCase)
            n_faults = 100;
            mf = MultiFault();
            mf = mf.initialize(n_faults);
            fault_length = length(mf.faults(1).y); 
            depth_dependent_friction = 0.6*ones(fault_length, 1) + rand(fault_length, 1)*0.1;
            % set the same depth-dependent input parameter for all faults
            mf.set_depth_dependent_input_parameter('f_s', repmat({depth_dependent_friction}, n_faults, 1));
            depth_dependent_friction = 0.5*ones(fault_length, 1);
            % set depth-dependent input parameter for one of the faults
            mf.set_depth_dependent_input_parameter('f_s', depth_dependent_friction, 2);
            f_s_for_check = mf.get_depth_dependent_input_parameter('f_s');
            testCase.verifyEqual(f_s_for_check{2}(1), 0.5);
            % different method to set depth-dependent input parameter
            mf.faults(2).set_depth_dependent_input_parameter('f_s', 0.55);
            f_s_for_check = mf.get_depth_dependent_input_parameter('f_s');
            testCase.verifyEqual(f_s_for_check{2}(1), 0.55);
        end

    end
end