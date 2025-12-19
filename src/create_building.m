function building = create_building(name, story_heights, wall_thickness, density, E, geom, nonlinear)
% CREATE_BUILDING Build a Building object from simple geometry and material data.
% The function:
%  - validates inputs
%  - removes floors with NaN thickness (treated as missing)
%  - expands scalar inputs to vectors
%  - computes per-floor mass and stiffness using get_m and get_k
%  - returns a Building instance
%
% All quantities are SI. Nonlinear must contain fields:
%   yield_drift, ultimate_drift, residual_strength, degradation_rate

    if nargin < 1 || isempty(name)
        name = "Building";
    end

    story_heights  = story_heights(:);
    wall_thickness = wall_thickness(:);

    if length(story_heights) ~= length(wall_thickness)
        error('story_heights and wall_thickness must have the same initial length.');
    end

    n_orig = length(story_heights);

    if isscalar(density), density = density * ones(n_orig,1); end
    if isscalar(E),       E       = E       * ones(n_orig,1); end
    density = density(:);
    E       = E(:);

    if any([length(density), length(E)] ~= n_orig)
        error('density and E must be scalar or vectors of same length as story_heights.');
    end

    % Remove floors with NaN thickness (treated as non-existent)
    valid = ~isnan(wall_thickness);
    if ~any(valid)
        error('All floors have NaN wall_thickness; nothing to build.');
    end

    story_heights  = story_heights(valid);
    wall_thickness = wall_thickness(valid);
    density        = density(valid);
    E              = E(valid);

    n = numel(story_heights); % updated number of floors

    % Validate geometry struct and expand scalar entries to vectors
    if ~isstruct(geom)
        error('geom must be a struct with fields wall_width, slab_thickness, internal_fraction (optional: floor_area, perimeter).');
    end
    if ~isfield(geom,'wall_width')
        error('geom must contain field ''wall_width''.');
    end
    if ~isfield(geom,'slab_thickness')
        error('geom must contain field ''slab_thickness''.');
    end
    if ~isfield(geom,'internal_fraction')
        error('geom must contain field ''internal_fraction''.');
    end

    if isscalar(geom.slab_thickness)
        slab_thickness = geom.slab_thickness * ones(n,1);
    else
        slab_thickness = geom.slab_thickness(:);
        if length(slab_thickness) ~= n
            error('geom.slab_thickness must be scalar or length n (after removing NaNs).');
        end
    end

    if isscalar(geom.wall_width)
        wall_width = geom.wall_width * ones(n,1);
    else
        wall_width = geom.wall_width(:);
        if length(wall_width) ~= n
            error('geom.wall_width must be scalar or length n (after removing NaNs).');
        end
    end

    if isfield(geom,'floor_area') && ~isempty(geom.floor_area)
        if isscalar(geom.floor_area)
            floor_area = geom.floor_area * ones(n,1);
        else
            floor_area = geom.floor_area(:);
            if length(floor_area) ~= n
                error('geom.floor_area must be scalar or length n.');
            end
        end
    else
        floor_area = wall_width .* wall_width; % assume square plan
    end

    if isfield(geom,'perimeter') && ~isempty(geom.perimeter)
        if isscalar(geom.perimeter)
            perimeter = geom.perimeter * ones(n,1);
        else
            perimeter = geom.perimeter(:);
            if length(perimeter) ~= n
                error('geom.perimeter must be scalar or length n.');
            end
        end
    else
        perimeter = 4 .* wall_width;
    end

    internal_fraction = geom.internal_fraction;
    if isscalar(internal_fraction)
        internal_fraction = internal_fraction * ones(n,1);
    else
        internal_fraction = internal_fraction(:);
        if length(internal_fraction) ~= n
            error('geom.internal_fraction must be scalar or length n.');
        end
    end

    % Compute per-floor mass and stiffness
    m = get_m(wall_thickness, density, story_heights, perimeter, floor_area, slab_thickness, internal_fraction);
    k = get_k(wall_thickness, E, story_heights, wall_width);

    if numel(m) ~= n || numel(k) ~= n
        error('Computed per-floor mass or stiffness length does not match number of valid stories.');
    end

    % Verify nonlinear struct has the required fields
    reqNL = {'yield_drift','ultimate_drift','residual_strength','degradation_rate'};
    for ii = 1:numel(reqNL)
        if ~isfield(nonlinear, reqNL{ii})
            error('nonlinear struct must contain field ''%s''.', reqNL{ii});
        end
    end

    building = Building(name, story_heights, m, k, nonlinear.yield_drift, nonlinear.ultimate_drift, nonlinear.residual_strength, nonlinear.degradation_rate);
end