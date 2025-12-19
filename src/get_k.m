function k = get_k(thickness, E, height, wall_width)
%GET_K Computes lateral stiffness for each floor
%
% INPUTS (all vectors of same length n_floors):
%   thickness   : wall thickness per floor [m]
%   E           : Young's modulus per floor [Pa]
%   height      : floor height per floor [m]
%   wall_width  : wall width per floor [m]
%
% OUTPUT:
%   k           : floor stiffness vector [N/m]

    % ensure column vectors
    thickness   = thickness(:);
    E           = E(:);
    height      = height(:);
    wall_width  = wall_width(:);

    % check consistency
    n = length(thickness);
    if any([length(E), length(height), length(wall_width)] ~= n)
        error('All input vectors must have the same length.');
    end

    % remove non-existing floors (NaN thickness)
    valid = ~isnan(thickness);

    thickness   = thickness(valid);
    E           = E(valid);
    height      = height(valid);
    wall_width  = wall_width(valid);

    % second moment of area
    I = wall_width .* thickness.^3 / 12;

    % lateral stiffness per floor
    k = 2 .* (12 .* E .* I ./ height.^3);
end
