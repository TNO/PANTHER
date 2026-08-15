classdef TestPantherParameterList < matlab.unittest.TestCase
    % Verify PantherParameterList creates independent parameter structs
    methods (Test)
        function uniqueParams(testCase)
            p1 = PantherParameterList();
            p2 = PantherParameterList();
            % Changing one shouldn't affect the other
            original = p2.young.value;
            p1.young.value = original + 1;
            testCase.verifyNotEqual(p1.young.value, p2.young.value);
        end
    end
end
