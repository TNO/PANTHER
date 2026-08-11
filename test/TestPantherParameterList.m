classdef TestPantherParameterList < matlab.unittest.TestCase
    % Verify PantherParameterList creates independent PantherParam instances
    methods (Test)
        function uniqueParams(testCase)
            p1 = PantherParameterList();
            p2 = PantherParameterList();
            % Ensure the parameter handles are distinct
            testCase.verifyFalse(p1.young == p2.young);
            % Changing one shouldn't affect the other
            original = p2.young.value;
            p1.young.value = original + 1;
            testCase.verifyNotEqual(p1.young.value, p2.young.value);
        end
    end
end
