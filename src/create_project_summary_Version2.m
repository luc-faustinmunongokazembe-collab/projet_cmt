function create_project_summary(projectRoot)
% CREATE_PROJECT_SUMMARY Create a compact tabular project summary for all buildings.
% Writes results/project_summary.txt with per-building rows:
% Country | Building | Stories | MaxDisp (mm) | StiffRed (%) | State
%
% The routine prefers explicit per-floor stiffness-reduction arrays saved
% in the individual MAT files; if those are not present it falls back to
% reasonable alternatives and remains silent on missing/partial data.

    if nargin < 1 || isempty(projectRoot)
        projectRoot = pwd();
    end

    dataDir = fullfile(projectRoot, 'data');
    resultsRoot = fullfile(projectRoot, 'results');

    if ~exist(resultsRoot, 'dir')
        error('Results folder not found: %s. Run simulations first.', resultsRoot);
    end

    % Try to find an AT2 file (best-effort, prefer earthquake_data.AT2)
    at2Path = '';
    if exist(dataDir,'dir')
        pref = 'earthquake_data.AT2';
        files = [dir(fullfile(dataDir,'*.AT2')); dir(fullfile(dataDir,'*.at2'))];
        if ~isempty(files)
            idx = find(strcmpi({files.name}, pref),1);
            if isempty(idx), idx = 1; end
            at2Path = fullfile(dataDir, files(idx).name);
        end
    end
    if isempty(at2Path)
        files = [dir('*.AT2'); dir('*.at2')];
        if ~isempty(files)
            at2Path = fullfile(pwd, files(1).name);
        end
    end

    at2_dt = NaN; at2_npts = NaN; at2_PGA_g = NaN;
    if ~isempty(at2Path)
        try
            [ag_g, dt, npts] = readAT2(at2Path);
            at2_dt = dt; at2_npts = npts; at2_PGA_g = max(abs(ag_g));
        catch
            % keep NaNs on failure, but do not produce console warnings
        end
    end

    outPath = fullfile(resultsRoot, 'project_summary.txt');
    fid = fopen(outPath, 'w');
    if fid == -1
        error('Could not open %s for writing', outPath);
    end

    % Header
    if ~isempty(at2Path) && ~isnan(at2_dt)
        fprintf(fid, 'AT2=%s, dt=%.6g, npts=%d, PGA_g=%.3f\n\n', at2Path, at2_dt, at2_npts, at2_PGA_g);
    elseif ~isempty(at2Path)
        fprintf(fid, 'AT2=%s\n\n', at2Path);
    else
        fprintf(fid, 'AT2=NONE\n\n');
    end

    % Table column widths
    wCountry = 15;
    wBuilding = 20;
    wStories = 6;
    wMaxDisp = 14;
    wStiff = 12;
    wState = 12;

    fprintf(fid, '%-*s %-*s %-*s %-*s %-*s %-*s\n', ...
        wCountry, 'Country', wBuilding, 'Building', wStories, 'Stories', ...
        wMaxDisp, 'MaxDisp (mm)', wStiff, 'StiffRed (%)', wState, 'State');
    totalWidth = wCountry + 1 + wBuilding + 1 + wStories + 1 + wMaxDisp + 1 + wStiff + 1 + wState;
    fprintf(fid, '%s\n', repmat('-', 1, totalWidth));

    % Iterate through results folders
    countryDirs = dir(resultsRoot);
    countryDirs = countryDirs([countryDirs.isdir] & ~ismember({countryDirs.name},{'.','..'}));
    for ci = 1:numel(countryDirs)
        countryName = countryDirs(ci).name;
        countryPath = fullfile(resultsRoot, countryName);
        buildingDirs = dir(countryPath);
        buildingDirs = buildingDirs([buildingDirs.isdir] & ~ismember({buildingDirs.name},{'.','..'}));

        for bi = 1:numel(buildingDirs)
            bname = buildingDirs(bi).name;
            bdir = fullfile(countryPath, bname);

            % Find newest mat file for the building
            matFiles = dir(fullfile(bdir, '*_mdof_results_building*.mat'));
            if isempty(matFiles)
                matFiles = dir(fullfile(bdir, '*.mat'));
            end
            matFull = '';
            if ~isempty(matFiles)
                [~, idxNewest] = max([matFiles.datenum]);
                matFull = fullfile(matFiles(idxNewest).folder, matFiles(idxNewest).name);
            end

            % Default metrics
            maxDisp_mm_overall = NaN;
            avg_stiff_red_pct = NaN;
            stateLabel = 'NoData';
            nStories = NaN;

            if ~isempty(matFull)
                try
                    S = load(matFull);

                    % Displacement metrics
                    if isfield(S,'u')
                        u = S.u;
                        nStories = size(u,1);
                        maxDisp_mm = max(abs(u),[],2) * 1000;
                        maxDisp_mm_overall = max(maxDisp_mm);
                    end

                    % Prefer explicit per-floor stiffness reductions
                    perFloorStiffRed_pct = [];
                    if isfield(S,'perFloorStiffRed_pct') && ~isempty(S.perFloorStiffRed_pct)
                        perFloorStiffRed_pct = S.perFloorStiffRed_pct(:);
                        avg_stiff_red_pct = mean(perFloorStiffRed_pct,'omitnan');
                    else
                        if isfield(S,'final_stiffness') && isfield(S,'k_elastic')
                            final_k = S.final_stiffness(:);
                            k_el = S.k_elastic(:);
                            if numel(final_k) == numel(k_el)
                                perFloorStiffRed_pct = (1 - final_k ./ k_el) * 100;
                                avg_stiff_red_pct = mean(perFloorStiffRed_pct,'omitnan');
                            end
                        else
                            if isfield(S,'final_stiffness') && isfield(S,'building')
                                final_k = S.final_stiffness(:);
                                k_el = [];
                                try
                                    bld = S.building;
                                    if isobject(bld) && isprop(bld,'k_elastic')
                                        k_el = bld.k_elastic(:);
                                    elseif isstruct(bld) && isfield(bld,'k_elastic')
                                        k_el = bld.k_elastic(:);
                                    end
                                catch
                                    k_el = [];
                                end
                                if ~isempty(k_el) && numel(final_k) == numel(k_el)
                                    perFloorStiffRed_pct = (1 - final_k ./ k_el) * 100;
                                    avg_stiff_red_pct = mean(perFloorStiffRed_pct,'omitnan');
                                end
                            end
                        end
                    end

                    % Determine compact state using available info
                    max_drift = NaN; yield_drift = NaN;
                    if isfield(S,'max_drift_per_story')
                        max_drift = max(S.max_drift_per_story(:));
                    elseif isfield(S,'story_drifts')
                        max_drift = max(max(abs(S.story_drifts)));
                    end
                    if isfield(S,'building')
                        try
                            bld = S.building;
                            if isobject(bld) && isprop(bld,'yield_drift')
                                yield_drift = bld.yield_drift;
                            elseif isstruct(bld) && isfield(bld,'yield_drift')
                                yield_drift = bld.yield_drift;
                            end
                        catch
                            yield_drift = NaN;
                        end
                    end

                    % State decision tree
                    if ~isempty(perFloorStiffRed_pct)
                        k_ratio = 1 - perFloorStiffRed_pct/100;
                        if any(k_ratio < 0.10)
                            stateLabel = 'Collapse';
                        elseif any(k_ratio < 0.30)
                            stateLabel = 'Severe';
                        elseif any(k_ratio < 0.60)
                            stateLabel = 'Moderate';
                        elseif ~isnan(max_drift) && ~isnan(yield_drift) && max_drift > yield_drift
                            stateLabel = 'Yielding';
                        else
                            stateLabel = 'Elastic';
                        end
                    elseif ~isnan(avg_stiff_red_pct)
                        if avg_stiff_red_pct > 70
                            stateLabel = 'Severe';
                        elseif avg_stiff_red_pct > 30
                            stateLabel = 'Moderate';
                        elseif ~isnan(max_drift) && ~isnan(yield_drift) && max_drift > yield_drift
                            stateLabel = 'Yielding';
                        else
                            stateLabel = 'Elastic';
                        end
                    elseif ~isnan(max_drift) && ~isnan(yield_drift) && max_drift > yield_drift
                        stateLabel = 'Yielding';
                    else
                        stateLabel = 'NoData';
                    end

                catch
                    % if parsing/loading fails, keep NoData silently
                    stateLabel = 'NoData';
                end
            else
                % try summary text as last resort, but stay silent
                txtFiles = dir(fullfile(bdir, '*_summary.txt'));
                if ~isempty(txtFiles)
                    try
                        txt = fileread(fullfile(bdir, txtFiles(1).name));
                        if contains(lower(txt), 'collapsed') || contains(lower(txt),'collapse')
                            stateLabel = 'Collapse';
                        elseif contains(lower(txt),'severe')
                            stateLabel = 'Severe';
                        elseif contains(lower(txt),'moderate')
                            stateLabel = 'Moderate';
                        elseif contains(lower(txt),'yield')
                            stateLabel = 'Yielding';
                        else
                            stateLabel = 'NoData';
                        end
                    catch
                        stateLabel = 'NoData';
                    end
                end
            end

            if isnan(maxDisp_mm_overall)
                maxDispStr = 'NA';
            else
                maxDispStr = sprintf('%.1f', maxDisp_mm_overall);
            end
            if isnan(avg_stiff_red_pct)
                stiffStr = 'NA';
            else
                stiffStr = sprintf('%.1f', avg_stiff_red_pct);
            end
            if isnan(nStories)
                nStoriesStr = 'NA';
            else
                nStoriesStr = sprintf('%d', nStories);
            end

            countrySafe = strrep(countryName, ',', ';');
            bnameSafe = strrep(bname, ',', ';');
            fprintf(fid, '%-*s %-*s %-*s %-*s %-*s %-*s\n', ...
                wCountry, trunc(countrySafe, wCountry), ...
                wBuilding, trunc(bnameSafe, wBuilding), ...
                wStories, nStoriesStr, ...
                wMaxDisp, maxDispStr, ...
                wStiff, stiffStr, ...
                wState, trunc(stateLabel, wState));
        end
    end

    fclose(fid);

    % Do not open the file automatically; it has been written to disk.

    function s = trunc(str, maxlen)
        if isempty(str)
            s = 'NA';
            return;
        end
        if length(str) <= maxlen
            s = str;
        else
            if maxlen > 3
                s = [str(1:maxlen-3) '...'];
            else
                s = str(1:maxlen);
            end
        end
    end
end