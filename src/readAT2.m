function [ag_g, dt, npts] = readAT2(filename)
%READAT2 Robust reader for PEER .AT2 files
%   [ag_g, dt, npts] = readAT2(filename)
%   - ag_g returned in units of 'g' (no conversion to m/s^2)
%   - dt returned in seconds
%   - npts returned is the number of samples actually returned (may differ from header)
%
% Behavior:
% - Searches first 100 lines for NPTS and DT tokens (robust to slight format variations).
% - Reads numeric tokens after header using textscan('%f'), trims or warns if counts mismatch.
% - If header NPTS is missing, npts is set to the number of numeric tokens read.

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    % Read up to first 100 lines to find NPTS and DT (don't assume fixed header length)
    maxHeaderLines = 100;
    headerLines = cell(maxHeaderLines,1);
    npts_header = [];
    dt_header = [];
    headerLineIdx = 0;
    for i = 1:maxHeaderLines
        tline = fgetl(fid);
        if ~ischar(tline)
            break; % EOF reached early
        end
        headerLines{i} = tline;
        % try to extract NPTS and DT tokens from this line
        npts_match = regexp(tline, 'NPTS\s*=\s*(\d+)', 'tokens', 'once');
        dt_match   = regexp(tline, 'DT\s*=\s*([-\d\.Ee+]+)', 'tokens', 'once');
        if ~isempty(npts_match)
            npts_header = str2double(npts_match{1});
        end
        if ~isempty(dt_match)
            dt_header = str2double(dt_match{1});
        end
        if ~isempty(npts_header) && ~isempty(dt_header)
            headerLineIdx = i;
            break; % found both, stop scanning header
        end
    end

    % If not found both, attempt to find them anywhere in headerLines already read
    if isempty(npts_header) || isempty(dt_header)
        for i = 1:maxHeaderLines
            if isempty(headerLines{i}), break; end
            tline = headerLines{i};
            if isempty(npts_header)
                npts_match = regexp(tline, 'NPTS\s*=\s*(\d+)', 'tokens', 'once');
                if ~isempty(npts_match), npts_header = str2double(npts_match{1}); end
            end
            if isempty(dt_header)
                dt_match = regexp(tline, 'DT\s*=\s*([-\d\.Ee+]+)', 'tokens', 'once');
                if ~isempty(dt_match), dt_header = str2double(dt_match{1}); end
            end
            if ~isempty(npts_header) && ~isempty(dt_header)
                headerLineIdx = i;
                break;
            end
        end
    end

    % If we still don't have header values, set headerLineIdx to last read line so we
    % start reading numeric tokens after those lines
    if headerLineIdx == 0
        % find last non-empty header line read
        lastNonEmpty = find(~cellfun(@isempty, headerLines), 1, 'last');
        if isempty(lastNonEmpty)
            lastNonEmpty = 0;
        end
        headerLineIdx = lastNonEmpty;
    end

    % Now position the file indicator to just after the header we scanned.
    % Re-open file and skip headerLineIdx lines to ensure consistent positioning.
    frewind(fid);
    for i = 1:headerLineIdx
        if ~ischar(fgetl(fid)), break; end
    end

    % Read numeric tokens from rest of file robustly
    dataCells = textscan(fid, '%f', 'CollectOutput', true);
    fclose(fid);

    if isempty(dataCells) || isempty(dataCells{1})
        error('No numeric acceleration data found in file: %s', filename);
    end

    ag_data = dataCells{1}(:); % column vector of numeric tokens

    % If header specified NPTS, trim or warn accordingly
    if ~isempty(npts_header)
        if length(ag_data) < npts_header
            warning('Header NPTS=%d but read %d numeric samples; using actual sample count.', npts_header, length(ag_data));
            npts = length(ag_data);
        elseif length(ag_data) > npts_header
            % More data than header says: trim to NPTS and warn (common if file contains trailing values)
            warning('Header NPTS=%d but read %d numeric samples; trimming to header count.', npts_header, length(ag_data));
            ag_data = ag_data(1:npts_header);
            npts = npts_header;
        else
            npts = npts_header;
        end
    else
        % No header NPTS found -> use actual count
        npts = length(ag_data);
        warning('No NPTS found in header; using actual sample count = %d.', npts);
    end

    % If header dt found use it, else try to extract dt from header last line or error
    if isempty(dt_header)
        % try last header line for two numbers (NPTS DT) e.g. fallback
        lastHeaderLine = headerLines{headerLineIdx};
        if ischar(lastHeaderLine)
            tmp = sscanf(lastHeaderLine, '%f %f', [1, Inf]);
            if length(tmp) >= 2
                dt_header = tmp(2);
            end
        end
    end
    if isempty(dt_header)
        error('Could not parse DT from header of %s. Please check file format.', filename);
    end
    dt = dt_header;

    % Return acceleration in units of g (do NOT convert here)
    ag_g = ag_data(:);
end