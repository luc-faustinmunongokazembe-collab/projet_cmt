classdef Building
    properties
        Name            % Name of the building
        story_heights   % Heights of each story [m]
        m               % Mass of each floor [kg]
        k_elastic       % Elastic stiffness of each floor [N/m]
        
        % Nonlinear progressive degradation parameters
        yield_drift     % Yield drift ratio
        ultimate_drift  % Ultimate drift ratio
        residual_strength % Residual strength (percentage of initial strength)
        degradation_rate  % Rate of stiffness degradation
    end
    
    methods
        function obj = Building(Name, story_heights, m, k_elastic, yield_drift, ultimate_drift, residual_strength, degradation_rate)
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