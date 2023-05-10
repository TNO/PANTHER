classdef GreensCalculator
    % GreensCalculator Combines Green function components for runs with diffusion 
    % runs , transforms to shear and normal stress. 

    properties
        a
        b
        Gshear
        Gnorm_FW
        Gshear_HW
        Gnorm_HW
        numerical_correction = 1e-9
    end

    methods
        function self = GreensCalculator(y, dx, params, pressure, temperature)
            

        end

        function self = get_components_uniform(self)
        % jjj
        end

        % steps for non-uniform
        % - recalculate a and b for each slice
        % - compute Greens function for each slice
        % - multiply slice Green function with C (=dP, elastic)
        % - sum all slices



    end

end