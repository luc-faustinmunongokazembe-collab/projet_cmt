function run_simulations()
% RUN_SIMULATIONS Run all country simulations and generate a project summary.
% Steps:
%  1) Locate project root (parent folder containing 'src')
%  2) Add src to the path if necessary
%  3) Attempt to compile the MEX solver if not already present
%  4) Ensure results/ is present
%  5) Load countries list and run each country's simulations
%  6) Create a consolidated project summary file
%
% The function is conservative about console output and avoids printing
% progress messages or non-fatal warnings.

    cur = pwd();
    projectRoot = '';
    while true
        if exist(fullfile(cur,'src'),'dir')
            projectRoot = cur;
            break;
        end
        parent = fileparts(cur);
        if isempty(parent) || strcmp(parent, cur)
            break;
        end
        cur = parent;
    end
    if isempty(projectRoot)
        error('Could not locate project root (folder containing ''src''). Please run from within the project tree.');
    end

    srcDir = fullfile(projectRoot, 'src');
    if ~exist(srcDir, 'dir')
        error('src directory not found at expected location: %s', srcDir);
    end
    if isempty(strfind(path, [srcDir pathsep])) && isempty(strfind(path, [pathsep srcDir]))
        addpath(srcDir);
    end

    mexName = 'mdof_degrade_mex';
    mexFullPath = fullfile(srcDir, [mexName, '.', mexext]);
    if ~exist(mexFullPath, 'file')
        oldpwd = pwd();
        try
            cd(srcDir);
            mex mdof_degrade_mex.c;
        catch
            % If compilation fails, continue quietly; user may compile manually.
        end
        cd(oldpwd);
    end

    resultsRoot = fullfile(projectRoot, 'results');
    if ~exist(resultsRoot, 'dir')
        mkdir(resultsRoot);
    end

    try
        countries_list = countries();
    catch ME
        error('Failed to load countries(): %s', ME.message);
    end

    for i = 1:length(countries_list)
        try
            countries_list(i).run_simulations();
        catch
            % continue to next country on failure without printing
        end
    end

    try
        create_project_summary(projectRoot);
    catch
        % ignore summary creation failures quietly
    end
end