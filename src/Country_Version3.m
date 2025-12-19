classdef Country
    % Country Holds a collection of Building objects and runs simulations.
    % The run_simulations method runs the MDOF solver for every building and
    % saves each building's outputs into results/<CountryName>/<BuildingName>/.
    properties
        Name       % Country name
        Buildings  % Array of Building objects
    end

    methods
        function obj = Country(name, buildings)
            % Create a Country with a name and a list of buildings.
            obj.Name = name;
            obj.Buildings = buildings;
        end

        function run_simulations(obj)
            % Run simulations for every building and save results files.
            %
            % The routine creates the country/result folders if needed and
            % calls the central driver run_mdof_mex for each building.
            country_dir = fullfile('results', obj.Name);
            if ~exist(country_dir, 'dir')
                mkdir(country_dir);
            end

            for i = 1:length(obj.Buildings)
                build = obj.Buildings(i);

                building_dir = fullfile(country_dir, build.Name);
                if ~exist(building_dir, 'dir')
                    mkdir(building_dir);
                end

                % Run the solver (the solver will save outputs into the building folder).
                run_mdof_mex(build, obj.Name, i);

                % Save a small reference copy of the Building object.
                save(fullfile(building_dir, [build.Name, '_results.mat']), 'build');
            end
        end
    end
end