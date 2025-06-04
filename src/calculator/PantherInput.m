function [instance_of_pantherinpu] = PantherInput()
    % temporary function to redirect to refactored class name
    % FaultSurfaceCalculator. Will be removed after v0.8.6
    instance_of_pantherinpu = PantherAnalysis();
    warning(['Class PantherInput has been refactored and will be removed after v0.8.6.',...
        'Use PantherAnalysis instead']);
end
