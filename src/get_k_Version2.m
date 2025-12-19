function k = get_k(thickness, E, height, wall_width)
% GET_K Compute lateral (shear) stiffness for each story.
% Inputs must be vectors of the same length (thickness, E, height, wall_width).
% Floors with NaN thickness are treated as non-existing and are removed.

    thickness   = thickness(:);
    E           = E(:);
    height      = height(:);
    wall_width  = wall_width(:);

    n = length(thickness);
    if any([length(E), length(height), length(wall_width)] ~= n)
        error('All input vectors must have the same length.');
    end

    valid = ~isnan(thickness);

    thickness   = thickness(valid);
    E           = E(valid);
    height      = height(valid);
    wall_width  = wall_width(valid);

    % Second moment of area for wall element (approximate)
    I = wall_width .* thickness.^3 / 12;

    % Shear-based lateral stiffness for each floor (two walls assumed)
    k = 2 .* (12 .* E .* I ./ height.^3);
end