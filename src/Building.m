classdef Building
    % Building Simple container for per-building parameters used by the solver.
    % This class stores geometric, mass/stiffness and nonlinear degradation
    % parameters for a multi-story building. It is intentionally lightweight:
    % the numerical routines live in separate functions and the compiled MEX.
    properties
        Name            % Name of the building
        story_heights   % Heights of each story [m]
        m               % Mass of each floor [kg]
        k_elastic       % Elastic stiffness of each floor [N/m]
        yield_drift     % Yield drift ratio (nonlinear parameter)
        ultimate_drift  % Ultimate drift ratio (nonlinear parameter)
        residual_strength % Residual strength (fraction of initial stiffness)
        degradation_rate  % Rate controlling stiffness degradation
    end

    methods
        function obj = Building(Name, story_heights, m, k_elastic, yield_drift, ultimate_drift, residual_strength, degradation_rate)
            % Construct a Building object.
            obj.Name = Name;
            obj.story_heights = story_heights;
            obj.m = m;
            obj.k_elastic = k_elastic;
            obj.yield_drift = yield_drift;
            obj.ultimate_drift = ultimate_drift;
            obj.residual_strength = residual_strength;
            obj.degradation_rate = degradation_rate;
        end
    end
end