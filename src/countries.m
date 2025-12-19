function countries = countries()
%COUNTRIES Defines all countries and their building inventories
%   Uses Option A: floor_area and perimeter are computed from wall_width

    % ---------- ANCIENT JAPAN (wood)----------
    geom_jp1.wall_width = 10;         % m per floor
    geom_jp1.slab_thickness = 0.1; % floor-specific slab thickness
    geom_jp1.internal_fraction = 0.15;          % fraction of walls + slab mass
    
    geom_jp2.wall_width = 6;         % m per floor
    geom_jp2.slab_thickness = geom_jp1.slab_thickness; % floor-specific slab thickness
    geom_jp2.internal_fraction = geom_jp1.internal_fraction;          % fraction of walls + slab mass
    
    geom_jp3.wall_width = 8;         % m per floor
    geom_jp3.slab_thickness = geom_jp1.slab_thickness; % floor-specific slab thickness
    geom_jp3.internal_fraction = geom_jp1.internal_fraction;          % fraction of walls + slab mass
    
    geom_jp4.wall_width = 2.5;         % m per floor
    geom_jp4.slab_thickness = geom_jp1.slab_thickness; % floor-specific slab thickness
    geom_jp4.internal_fraction = geom_jp1.internal_fraction;          % fraction of walls + slab mass
    
    thickness_jp = .15;

    nl_jp.yield_drift = 0.005;
    nl_jp.ultimate_drift = 0.06;
    nl_jp.residual_strength = 0.4;
    nl_jp.degradation_rate = 0.6;

    E_jp = 10.5e9;
    rho_jp = 500;

    building_jp1 = create_building('JP_Building_1',[3;3;2.8;2.8], thickness_jp * ones(4,1), rho_jp, E_jp, geom_jp1, nl_jp);

    building_jp2 = create_building('JP_Building_2',[3;3], thickness_jp * ones(2,1), rho_jp, E_jp, geom_jp2, nl_jp);

    building_jp3 = create_building('JP_Building_3',[3;2.8;2.8], thickness_jp * ones(3,1), rho_jp, E_jp, geom_jp3, nl_jp);

    building_jp4 = create_building('JP_Building_4',[1.7;1.5], thickness_jp * ones(2,1), rho_jp, E_jp, geom_jp4, nl_jp);

    japan = Country('Ancient Japan', [building_jp1, building_jp2,building_jp3,building_jp4]);

    % ---------- ANCIENT GREECE (stone)----------
    geom_gr1.wall_width = 8;                 % m per floor
    geom_gr1.slab_thickness = 0.2;       % floor-specific
    geom_gr1.internal_fraction = 0.05;

    geom_gr2.wall_width = 9;                 % m per floor
    geom_gr2.slab_thickness = 0.2;       % floor-specific
    geom_gr2.internal_fraction = 0.05;

    geom_gr3.wall_width = 15;                 % m per floor
    geom_gr3.slab_thickness = 0.2;       % floor-specific
    geom_gr3.internal_fraction = 0.05;

    geom_gr4.wall_width = 7;                 % m per floor
    geom_gr4.slab_thickness = 0.2;       % floor-specific
    geom_gr4.internal_fraction = 0.05;
    
    thickness_gr1 = .6;
    thickness_gr3 = .8;
    thickness_gr4 = .7;

    nl_gr.yield_drift = 0.0008;
    nl_gr.ultimate_drift = 0.004;
    nl_gr.residual_strength = 0.02;
    nl_gr.degradation_rate = 2.0;

    E_gr = 2.3e10;
    rho_gr = 2800;

    building_gr1 = create_building('GR_Building_1',[4;4], thickness_gr1 * ones(2,1), rho_gr, E_gr, geom_gr1, nl_gr);
    building_gr2 = create_building('GR_Building_2',[4;3.5;3], thickness_gr1 * ones(3,1), rho_gr, E_gr, geom_gr2, nl_gr);
    building_gr3 = create_building('GR_Building_3',5*ones(1,1), thickness_gr3*ones(1,1), rho_gr, E_gr, geom_gr3, nl_gr);
    building_gr4 = create_building('GR_Building_4',[3.5;3.5;3.5;3.5], thickness_gr4 * ones(4,1), rho_gr, E_gr, geom_gr4, nl_gr);

    greece = Country('Ancient Greece', [building_gr1, building_gr2, building_gr3,building_gr4]);
    
    % ---------- USA (Reinforced Concrete) ----------
    geom_us1.wall_width = 15;                 % m per floor
    geom_us1.slab_thickness = 0.18;       % floor-specific
    geom_us1.internal_fraction = 0.3;

    geom_us2.wall_width = 15;                 % m per floor
    geom_us2.slab_thickness = 0.18;       % floor-specific
    geom_us2.internal_fraction = 0.3;

    geom_us3.wall_width = 18;                 % m per floor
    geom_us3.slab_thickness = 0.18;       % floor-specific
    geom_us3.internal_fraction = 0.3;

    geom_us4.wall_width = 12;                 % m per floor
    geom_us4.slab_thickness = 0.18;       % floor-specific
    geom_us4.internal_fraction = 0.3;
    
    thickness_us = .2;

    nl_us.yield_drift = 0.004;
    nl_us.ultimate_drift = 0.04;
    nl_us.residual_strength = 0.3;
    nl_us.degradation_rate = 0.8;

    E_us = 30e9;
    rho_us = 2500;

    building_us1 = create_building('US_Building_1',[4.0; 3.8; 3.8; 3.5; 3.5], thickness_us * ones(5,1), rho_us, E_us, geom_us1, nl_us);
    building_us2 = create_building('US_Building_2',[4.0; 3.8; 3.5; 3.5], thickness_us * ones(4,1), rho_us, E_us, geom_us2, nl_us);
    building_us3 = create_building('US_Building_3',[4;3.5;3.5;3.5;3.5;3.5;3.5;3.5], thickness_us*ones(8,1), rho_us, E_us, geom_us3, nl_us);
    building_us4 = create_building('US_Building_4',[4;4], thickness_us * ones(2,1), rho_us, E_us, geom_us4, nl_us);

    usa = Country('USA', [building_us1, building_us2, building_us3,building_us4]);



    % ---------- Output list ----------
    countries = [japan, greece,usa];
end
