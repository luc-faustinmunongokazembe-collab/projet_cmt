function m = get_m(thickness, density, floor_height, effective_perimeter, floor_area, slab_thickness, internal_structure_fraction)
%GET_M Computes the mass of each floor of a building
%
% INPUTS (vectors must have the same length n_floors):
%   thickness                   : wall thickness per floor [m]
%   density                     : material density per floor [kg/m^3]
%   floor_height                : floor height per floor [m]
%   effective_perimeter         : effective wall perimeter per floor [m]
%   floor_area                  : floor area per floor [m^2]
%   slab_thickness              : slab thickness per floor [m]
%   internal_structure_fraction : fraction of (walls + slab) mass [-]
%
% OUTPUT:
%   m                           : floor mass vector [kg]

    % ensure column vectors
    thickness           = thickness(:);
    density             = density(:);
    floor_height        = floor_height(:);
    effective_perimeter = effective_perimeter(:);
    floor_area          = floor_area(:);
    slab_thickness      = slab_thickness(:);

    % check consistency
    n = length(thickness);
    if any([ length(density), length(floor_height),length(effective_perimeter), length(floor_area), length(slab_thickness) ] ~= n)
        error('All input vectors must have the same length.');
    end

    % remove non-existing floors
    valid = ~isnan(thickness);

    thickness           = thickness(valid);
    density             = density(valid);
    floor_height        = floor_height(valid);
    effective_perimeter = effective_perimeter(valid);
    floor_area          = floor_area(valid);
    slab_thickness      = slab_thickness(valid);

    % wall volume and mass
    wall_volume = effective_perimeter .* floor_height .* thickness;
    wall_mass   = wall_volume .* density;

    % slab volume and mass
    slab_volume = floor_area .* slab_thickness;
    slab_mass   = slab_volume .* density;

    % internal structural mass
    internal_mass = internal_structure_fraction .* (wall_mass + slab_mass);

    % total floor mass
    m = wall_mass + slab_mass + internal_mass;
end
