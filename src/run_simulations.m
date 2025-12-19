function run_simulations()
%RUN_SIMULATIONS Runs simulations for all countries and generates a project summary
% Locates project root (parent containing 'src'), ensures src on path,
% compiles MEX if needed, runs simulations, then creates and opens project_summary.txt.

    % 1) Locate project root (folder that contains 'src')
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

    % 2) Ensure src on path
    srcDir = fullfile(projectRoot, 'src');
    if ~exist(srcDir, 'dir')
        error('src directory not found at expected location: %s', srcDir);
    end
    if isempty(strfind(path, [srcDir pathsep])) && isempty(strfind(path, [pathsep srcDir]))
        addpath(srcDir);
    end

    % 3) Compile MEX if needed (in src so compiled binary sits next to .c)
    mexName = 'mdof_degrade_mex';
    mexFullPath = fullfile(srcDir, [mexName, '.', mexext]);
    if ~exist(mexFullPath, 'file')
        oldpwd = pwd();
        try
            cd(srcDir);
            fprintf('Compiling %s in %s ...\n', mexName, srcDir);
            mex mdof_degrade_mex.c;
            fprintf('Compilation finished.\n');
        catch ME
            warning('Could not compile %s: %s\nProceeding; make sure the compiled MEX is available in src.', mexName, ME.message);
        end
        cd(oldpwd);
    else
        fprintf('Found compiled MEX: %s\n', mexFullPath);
    end

    % 4) Ensure results folder exists at project root
    resultsRoot = fullfile(projectRoot, 'results');
    if ~exist(resultsRoot, 'dir')
        mkdir(resultsRoot);
    end

    % 5) Run simulations
    try
        countries_list = countries(); % function must be on path (src)
    catch ME
        error('Failed to load countries(): %s', ME.message);
    end

    for i = 1:length(countries_list)
        try
            fprintf('\n=== Running simulations for %s ===\n', countries_list(i).Name);
            countries_list(i).run_simulations();
        catch ME
            warning('Simulations for %s failed: %s', countries_list(i).Name, ME.message);
        end
    end

    % 6) Create consolidated project summary and open it in the MATLAB editor
    try
        create_project_summary(projectRoot);
    catch ME
        warning('Failed to create project summary: %s', ME.message);
    end
end