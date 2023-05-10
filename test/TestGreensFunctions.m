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
            yeval = get_fault_y(1, 500);
            xeval = yeval/(tan(dip*pi/180)) + dx;
            ycf = -300;
            xcf = ycf/(tan(dip*pi/180));
            greens_f = GreensFunctions(yeval);
            greens_f = greens_f.green_FW(xeval, yeval, dip, testCase.thick, throw, width, 0, 0 );
            mid_i = floor(length(yeval)/2);
            actual = greens_f.Gnorm_FW(mid_i);
            expected = -pi;
            testCase.verifyEqual(actual, expected , "RelTol", 0.01);     
        end
    end
end

%https://github.com/marketplace/actions/run-matlab-tests