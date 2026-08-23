classdef TestFaultParameterList < matlab.unittest.TestCase
    % Verify FaultParameterList creates independent parameter structs
    methods (Test)
        function uniqueParams(testCase)
            p1 = FaultParameterList();
            p2 = FaultParameterList();
            % Changing one shouldn't affect the other
            original = p2.young.value;
            p1.young.value = original + 1;
            testCase.verifyNotEqual(p1.young.value, p2.young.value);
        end
    end
end
