function m = get_m(thickness, density, floor_height, effective_perimeter, floor_area, slab_thickness, internal_structure_fraction)
% GET_M Compute the mass of each floor.
% The mass includes wall mass, slab mass and a fraction for internal structural mass.
% Floors with NaN thickness are treated as missing and removed.

    thickness           = thickness(:);
    density             = density(:);
    floor_height        = floor_height(:);
    effective_perimeter = effective_perimeter(:);
    floor_area          = floor_area(:);
    slab_thickness      = slab_thickness(:);

    n = length(thickness);
    if any([ length(density), length(floor_height), length(effective_perimeter), length(floor_area), length(slab_thickness) ] ~= n)
        error('All input vectors must have the same length.');
    end

    valid = ~isnan(thickness);

    thickness           = thickness(valid);
    density             = density(valid);
    floor_height        = floor_height(valid);
    effective_perimeter = effective_perimeter(valid);
    floor_area          = floor_area(valid);
    slab_thickness      = slab_thickness(valid);

    wall_volume = effective_perimeter .* floor_height .* thickness;
    wall_mass   = wall_volume .* density;

    slab_volume = floor_area .* slab_thickness;
    slab_mass   = slab_volume .* density;

    internal_mass = internal_structure_fraction .* (wall_mass + slab_mass);

    m = wall_mass + slab_mass + internal_mass;
end