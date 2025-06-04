function [instance_of_faultsurfacecalculator] = MultiFaultCalculator(n_pillars)
    % temporary function to redirect to refactored class name
    % FaultSurfaceCalculator. Will be removed after v0.8.6
    instance_of_faultsurfacecalculator = FaultSurfaceCalculator(n_pillars);
    warning(['Class MultiFaultCalculator has been refactored and will be removed after v0.8.6.',...
        'Use FaultSurfaceCalculator instead']);
end
