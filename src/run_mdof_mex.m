function run_mdof_mex(building, country, buildingIndex, showPlots)
% RUN_MDOF_MEX Run the compiled MEX solver for a building and write outputs.
% The function:
%  - discovers the ground motion file
%  - assembles M, K and Rayleigh damping
%  - calls mdof_degrade_mex
%  - post-processes results and writes MAT and text summaries and PNGs
%
% Non-essential console output and warnings have been removed so this
% function operates quietly (errors remain visible).

    if nargin < 4
        showPlots = false;
    end

    % Find project root (parent containing 'src') or default to cwd
    thisFile = mfilename('fullpath');
    if isempty(thisFile)
        thisDir = pwd();
    else
        thisDir = fileparts(thisFile);
    end

    [~, thisDirName] = fileparts(thisDir);
    if strcmpi(thisDirName, 'src')
        projectRoot = fileparts(thisDir);
    else
        projectRoot = '';
        cur = thisDir;
        while true
            if exist(fullfile(cur, 'src'), 'dir')
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
            parent = fileparts(thisDir);
            if exist(fullfile(parent, 'data'), 'dir') || exist(fullfile(parent, 'results'), 'dir')
                projectRoot = parent;
            else
                projectRoot = pwd();
            end
        end
    end

    dataDir = fullfile(projectRoot, 'data');
    resultsRoot = fullfile(projectRoot, 'results');
    if ~exist(resultsRoot, 'dir'), mkdir(resultsRoot); end

    story_heights = building.story_heights(:);
    n = length(story_heights);

    m = building.m(:);
    k_elastic = building.k_elastic(:);

    yield_drift       = building.yield_drift;
    ultimate_drift    = building.ultimate_drift;
    residual_strength = building.residual_strength;
    degradation_rate  = building.degradation_rate;

    % Discover ground motion file quietly
    preferredName = 'earthquake_data.AT2';
    candidates = [];
    if exist(dataDir, 'dir')
        candidates = [dir(fullfile(dataDir, '*.AT2')); dir(fullfile(dataDir, '*.at2'))];
    end
    if isempty(candidates)
        candidates = [dir('*.AT2'); dir('*.at2')];
        if ~isempty(candidates)
            dataDir = pwd();
        end
    end
    if isempty(candidates)
        error('Could not find any .AT2 files in %s or the current folder. Place your AT2 file in the data/ folder.', fullfile(projectRoot, 'data'));
    end

    chosen = [];
    for ii = 1:numel(candidates)
        if strcmpi(candidates(ii).name, preferredName)
            chosen = candidates(ii);
            break;
        end
    end
    if isempty(chosen)
        chosen = candidates(1);
    end

    filepath = fullfile(dataDir, chosen.name);

    % Assemble M and initial K, then compute Rayleigh damping
    M = diag(m);
    K_initial = build_shear_stiffness(k_elastic, n);

    lambda_all = sort(eig(K_initial, M));
    xi = 0.05;
    if length(lambda_all) >= 2 && all(lambda_all > 0)
        omega1 = sqrt(lambda_all(1));
        omega2 = sqrt(lambda_all(2));
        a0_ray = 2 * xi * omega1 * omega2 / (omega1 + omega2);
        a1_ray = 2 * xi / (omega1 + omega2);
        C = a0_ray * M + a1_ray * K_initial;
    else
        omega1 = sqrt(max(lambda_all(1), eps));
        C = 2 * xi * omega1 * M;
    end

    beta  = 1/4;
    gamma = 1/2;

    % Read ground motion (ag_g in units of g)
    [ag_g, dt, npts] = readAT2(filepath);
    ag = ag_g * 9.81;
    t  = (0:npts-1) * dt;
    PGA_g = max(abs(ag_g));

    % Prepare results directories
    results_dir = fullfile(resultsRoot, country);
    if ~exist(results_dir, 'dir'), mkdir(results_dir); end
    building_dir = fullfile(results_dir, building.Name);
    if ~exist(building_dir, 'dir'), mkdir(building_dir); end

    % Ensure compiled MEX is available
    if ~exist('mdof_degrade_mex', 'file')
        error('MEX function mdof_degrade_mex not found on path. Compile it (mex mdof_degrade_mex.c) and ensure src is on the MATLAB path.');
    end

    % Call solver
    [u, ud, udd, story_drifts, stiffness_history, deg_hist, plastic_hist, max_drift_ratio, yielded] = ...
        mdof_degrade_mex(M, C, k_elastic, story_heights, ag, dt, beta, gamma, ...
        yield_drift, ultimate_drift, residual_strength, degradation_rate);

    % Post-process shapes and compute per-floor stiffness reductions
    k_elastic = k_elastic(:);
    final_stiffness = stiffness_history(:, end);
    final_stiffness = final_stiffness(:);

    if numel(final_stiffness) ~= numel(k_elastic)
        if isscalar(final_stiffness)
            final_stiffness = repmat(final_stiffness, numel(k_elastic), 1);
        else
            final_stiffness = final_stiffness(:);
        end
    end

    perFloorStiffRed_pct = (k_elastic - final_stiffness) ./ k_elastic * 100;
    perFloorStiffRed_pct(~isfinite(perFloorStiffRed_pct)) = NaN;

    max_drift_per_story = max(max_drift_ratio, [], 2);
    stiffness_reduction = perFloorStiffRed_pct;
    plastic_drifts      = plastic_hist(:, end);
    colors = lines(max(1,n));

    % Save MAT
    matFile = fullfile(building_dir, sprintf('%s_mdof_results_building%d.mat', building.Name, buildingIndex));
    save(matFile, 'u','ud','udd','story_drifts','stiffness_history','deg_hist','plastic_hist','max_drift_ratio','yielded', ...
        't','dt','ag','ag_g','PGA_g','building','final_stiffness','k_elastic','perFloorStiffRed_pct','stiffness_reduction','max_drift_per_story','plastic_drifts');

    % Save displacement pages (failures handled silently)
    try
        perPage = 12;
        pages = ceil(n / perPage);
        for p = 1:pages
            first = (p-1)*perPage + 1;
            last = min(n, p*perPage);
            nplot = last - first + 1;
            fig = figure('Visible', ternary(showPlots,'on','off'), 'Position', [100 100 900 1200]);
            for s = first:last
                axIdx = s - first + 1;
                subplot(nplot, 1, axIdx);
                plot(t, u(s,:) * 1000, 'b-', 'LineWidth', 1.0); grid on;
                xlabel('Time (s)');
                ylabel('Disp (mm)');
                title(sprintf('Story %d', s));
                [~, idx_peak] = max(abs(u(s,:)));
                hold on;
                plot(t(idx_peak), u(s,idx_peak)*1000, 'ro', 'MarkerFaceColor', 'r');
                plot(t(end), u(s,end)*1000, 'k.', 'MarkerSize', 6);
                hold off;
            end
            sgtitle(sprintf('%s - Building %d (%s) - Floor Displacements (page %d/%d)', building.Name, buildingIndex, country, p, pages),'FontSize',14,'FontWeight','bold');
            outFile = fullfile(building_dir, sprintf('%s_figure_disp_page%d.png', building.Name, p));
            try
                exportgraphics(fig, outFile, 'Resolution', 300);
            catch
                try saveas(fig, outFile); catch, end
            end
            close(fig);
        end
    catch
        % continue quietly on plot generation errors
    end

    % Degradation plot (quiet on failures)
    try
        fig = figure('Visible', ternary(showPlots,'on','off'), 'Position', [200 200 1000 600]);
        hold on; grid on;
        for s = 1:n
            plot(t, deg_hist(s,:), 'LineWidth', 1.2, 'Color', colors(mod(s-1,size(colors,1))+1,:));
        end
        xlabel('Time (s)'); ylabel('Degradation index k(t)/k_initial');
        title(sprintf('%s - Building %d (%s) - Damage Evolution', building.Name, buildingIndex, country));
        legend(arrayfun(@(x) sprintf('Story %d', x), 1:n, 'UniformOutput', false), 'Location', 'bestoutside');
        outFile = fullfile(building_dir, sprintf('%s_degradation_index.png', building.Name));
        try exportgraphics(fig, outFile, 'Resolution', 300); catch, end
        close(fig);
    catch
        % continue quietly
    end

    % Write textual summary
    try
        summary_filename = fullfile(building_dir, sprintf('%s_building%d_summary.txt', building.Name, buildingIndex));
        fid = fopen(summary_filename,'w');
        if fid == -1, error('Could not open summary file for writing: %s', summary_filename); end

        fprintf(fid, 'EARTHQUAKE RECORD:\n');
        fprintf(fid, 'File: %s\n', filepath);
        fprintf(fid, 'PGA: %.3f g\n', PGA_g);
        fprintf(fid, 'Duration: %.1f s\n', t(end));
        fprintf(fid, 'Time step: %.4f s | Number of steps: %d\n\n', dt, length(t));

        fprintf(fid, 'DISPLACEMENT RESULTS:\n');
        fprintf(fid, '%-8s %-15s %-15s %-15s\n', 'Floor', 'Max Disp (mm)', 'Residual Disp (mm)', 'Drift Ratio (%)');
        fprintf(fid, '%s\n', repmat('-',1,60));
        for s = 1:n
            max_disp = max(abs(u(s,:))) * 1000;
            res_disp = u(s,end) * 1000;
            drift_ratio = max(abs(story_drifts(s,:))) * 100;
            fprintf(fid, '%-8d %-15.1f %-15.1f %-15.3f\n', s, max_disp, res_disp, drift_ratio);
        end

        fprintf(fid, '\nFINAL RESULTS:\n');
        fprintf(fid, '%-8s %-12s %-15s %-15s %-20s %-10s\n', 'Story', 'Max Drift(%)', 'Plastic Drift(%)', 'Stiffness(%)', 'Damage State', 'Yielded?');
        damage_labels = {'No Damage', 'Yielding', 'Moderate', 'Severe', 'Collapse'};
        final_damage_state = zeros(n,1);
        for s = 1:n
            k_ratio = stiffness_history(s,end)/k_elastic(s);
            max_d   = max_drift_per_story(s);
            if k_ratio < 0.10
                final_damage_state(s) = 4;
            elseif k_ratio < 0.30
                final_damage_state(s) = 3;
            elseif k_ratio < 0.60
                final_damage_state(s) = 2;
            elseif max_d > yield_drift
                final_damage_state(s) = 1;
            else
                final_damage_state(s) = 0;
            end
            fprintf(fid, '%-8d %-12.3f %-15.3f %-15.1f %-20s %-10s\n', s, ...
                max_d*100, plastic_hist(s,end)*100, k_ratio*100, damage_labels{final_damage_state(s)+1}, ternary(yielded(s),'YES','NO'));
        end

        fprintf(fid, '\nOVERALL BUILDING CONDITION:\n');
        if any(final_damage_state == 4)
            fprintf(fid, 'BUILDING COLLAPSED - Complete stiffness loss in one or more stories\n');
        elseif any(final_damage_state == 3)
            fprintf(fid, 'SEVERE DAMAGE - Significant stiffness reduction (>70%% loss)\n');
        elseif any(final_damage_state == 2)
            fprintf(fid, 'MODERATE DAMAGE - Noticeable stiffness reduction\n');
        elseif any(yielded)
            fprintf(fid, 'MINOR DAMAGE - Some yielding but stiffness largely intact\n');
        else
            fprintf(fid, 'BUILDING ELASTIC - No yielding, structure remains linear\n');
        end

        fprintf(fid, '\nNONLINEAR ANALYSIS STATISTICS:\n');
        fprintf(fid, 'Yielded stories: %d/%d (%.0f%%)\n', sum(yielded), n, sum(yielded)/n*100);
        fprintf(fid, 'Average stiffness reduction: %.1f%%\n', mean(stiffness_reduction,'omitnan'));
        fprintf(fid, 'Maximum plastic drift: %.3f%%\n', max(plastic_drifts)*100);

        fclose(fid);
    catch
        if exist('fid','var') && fid > 0, fclose(fid); end
    end
end

function K = build_shear_stiffness(k, n)
    K = zeros(n,n);
    for i = 1:n
        K(i,i) = K(i,i) + k(i);
        if i>1
            K(i,i-1) = -k(i);
            K(i-1,i) = -k(i);
            K(i-1,i-1) = K(i-1,i-1) + k(i);
        end
    end
end

function result = ternary(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end