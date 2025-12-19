function [ag_g, dt, npts] = readAT2(filename)
% READAT2 Read PEER .AT2 acceleration files (best-effort, robust parser).
% Returns acceleration in units of g, the time step dt (s) and number of points.
% The function is tolerant of some header variations. It does not print to console.

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    maxHeaderLines = 100;
    headerLines = cell(maxHeaderLines,1);
    npts_header = [];
    dt_header = [];
    headerLineIdx = 0;
    for i = 1:maxHeaderLines
        tline = fgetl(fid);
        if ~ischar(tline)
            break;
        end
        headerLines{i} = tline;
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
            break;
        end
    end

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

    if headerLineIdx == 0
        lastNonEmpty = find(~cellfun(@isempty, headerLines), 1, 'last');
        if isempty(lastNonEmpty)
            lastNonEmpty = 0;
        end
        headerLineIdx = lastNonEmpty;
    end

    frewind(fid);
    for i = 1:headerLineIdx
        if ~ischar(fgetl(fid)), break; end
    end

    dataCells = textscan(fid, '%f', 'CollectOutput', true);
    fclose(fid);

    if isempty(dataCells) || isempty(dataCells{1})
        error('No numeric acceleration data found in file: %s', filename);
    end

    ag_data = dataCells{1}(:);

    if ~isempty(npts_header)
        if length(ag_data) < npts_header
            npts = length(ag_data);
        elseif length(ag_data) > npts_header
            ag_data = ag_data(1:npts_header);
            npts = npts_header;
        else
            npts = npts_header;
        end
    else
        npts = length(ag_data);
    end

    if isempty(dt_header)
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

    ag_g = ag_data(:);
end