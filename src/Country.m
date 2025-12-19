classdef Country
    properties
        Name       % Country name
        Buildings  % Array of building objects
    end
    
    methods
        function obj = Country(name, buildings)
            obj.Name = name;
            obj.Buildings = buildings;  % Buildings will be an array of building objects
        end
        
        function run_simulations(obj)
            % Create a folder for the country if it doesn't exist
            country_dir = fullfile('results', obj.Name);
            if ~exist(country_dir, 'dir')
                mkdir(country_dir);
            end
            
            % Loop through each building in the country
            for i = 1:length(obj.Buildings)
                build = obj.Buildings(i);
                disp(['Running simulation for ', build.Name, ' in ', obj.Name]);
    
                % Create directory for building if it doesn't exist
                building_dir = fullfile(country_dir, build.Name);
                if ~exist(building_dir, 'dir')
                    mkdir(building_dir);
                end
    
                % Run the MDOF simulation for this building
                run_mdof_mex(build, obj.Name, i);
    
                % Save results (can be more specific if you only want to save results)
                save(fullfile(building_dir, [build.Name, '_results.mat']), 'build');
            end
        end
    end
end