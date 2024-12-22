classdef TestGreensFunctions < matlab.unittest.TestCase
    % TestGreensFunctions Unit test for the implementation of the 
    % Green's functions from Jansen et al. 2019
    % TODO: xx, yy, xy components
    % single side, addition. subtraction of 2 pi. offset along fault
    
    properties
        thick = 200;
    end
    
    methods (Test)
        function Gnorm_footwall(testCase)
            dip = 90;
            throw = 50;
            width = 1e9;
            dx = 1e-9;
            yeval = 0;
            xeval = yeval/(tan(dip*pi/180)) + dx;
            ycf = -300;
            xcf = ycf/(tan(dip*pi/180));
            greens_f = GreensFunctions(yeval);
            greens_f = greens_f.green_FW(xeval, yeval, dip, testCase.thick, throw, width, 0, 0 );
            actual = greens_f.Gnorm_FW(1);
            expected = -pi;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);     
        end
    end
end

%https://github.com/marketplace/actions/run-matlab-tests