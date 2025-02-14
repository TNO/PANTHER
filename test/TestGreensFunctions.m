classdef TestGreensFunctions < matlab.unittest.TestCase
    % TestGreensFunctions Unit test for the implementation of the 
    % Green's functions from Jansen et al. 2019
    
    properties
        thick = 200;
    end
    
    methods (Test)
        function Gnorm_footwall(testCase)
            % test for the footwall side, throw 50 m, thickness 200 m, 90
            % degree dip. Evaluation point within the reservoir
            dip = 90;
            throw = 50;
            width = 1e9;
            dx = 1e-9;
            yeval = 0;
            xeval = yeval/(tan(dip*pi/180)) + dx;
            greens_f = GreensFunctions(yeval);
            greens_f = greens_f.green_FW(xeval, yeval, dip, testCase.thick, throw, width, 0, 0 );
            actual = greens_f.Gnorm_FW(1);
            expected = -pi;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);     
        end
    end
end