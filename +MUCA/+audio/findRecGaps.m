function [dtGaps, dtRecs, badFiles] = findRecGaps(wavDir, varargin)
%
% Returns intervals that mark the occurrence of gaps between a series of
% audio files. May also plot the gaps if specified.
%
% SYNTAX:
%   [dtGaps, dtRecs, badFiles] = findRecGaps(wavDir)
%   [__] = findRecGaps(wavDir, recursive)
%   [__] = findRecGaps(__, 'Plot',value)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "wavDir" - path to folder containing WAV files
%   .......................................................................
%   "recursive" [OPTIONAL] - true/false value indicating if subfolders
%       should be searched or not. Default is false.
%   .......................................................................
%   "Tol" [OPTIONAL, NAME-VALUE PAIR] - number indicating the maximum
%       amount of time (in seconds) that can occur between successive audio
%       files before the difference is counted as a gap. Default is 0.
%   .......................................................................
%   "Plot" [OPTIONAL, NAME-VALUE PAIR] - true/false value indicating
%       whether or not to draw a plot to visualize the occurrence of gaps.
%       Default is false.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "dtGaps" - nGaps-by-2 matrix of datetime objects representing the start
%       and end times of gaps between recordings
%   .......................................................................
%   "dtRecs" - nRecs-by-2 matrix of datetime objects representing the start
%       and end times of each audio file in the time series
%   .......................................................................
%   "badFiles" - table with information on any files that failed to process
%   .......................................................................
%
% OUTPUT FILES:
%   <none>
%
% DEPENDENCIES:
%   MUCA.filepaths.listFiles
%   MUCA.audio.getRecTimes
%
%
% Written by Wilfried Beslin
% Last updated 2023-11-30, using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import MUCA.filepaths.listFiles
    import MUCA.audio.getRecTimes

    % input validation
    p = inputParser;
    
    % wavDir
    addRequired(p, 'wavDir', @isfolder)
    % recursive [OPTIONAL]
    defaultRecursive = false;
    validRecursive = @(arg) validateattributes(arg,{'logical'},{'scalar'});
    addOptional(p, 'recursive', defaultRecursive, validRecursive)
    % Tol [PARAMETER]
    defaultTol = 0;
    validTol = @(arg) validateattributes(arg,{'numeric'},{'scalar','nonnegative'});
    addParameter(p, 'Tol', defaultTol, validTol);
    % Plot [PARAMETER]
    defaultPlot = false;
    validPlot = @(arg) validateattributes(arg,{'logical'},{'scalar'});
    addParameter(p, 'Plot', defaultPlot, validPlot)
    
    parse(p,wavDir, varargin{:})
    recursive = p.Results.recursive;
    tol = seconds(p.Results.Tol);
    doPlot = p.Results.Plot;
    % end input parsing

    % get list of WAV files
    [wavPaths, wavNames] = listFiles(wavDir, 'wav', 'Recursive',recursive);
    numRecs = numel(wavNames);
    
    % get start and end times of each recording
    [dtStarts, dtStops] = getRecTimes(wavPaths);
    
    % identify problem files
    fileErrors = repmat({''}, numRecs, 1);
    noStartTimeStr = 'Unable to read file timestamp';
    noStopTimeStr = 'Could not determine file end time';
    noStartTime = isnat(dtStarts);
    noStopTime = isnat(dtStops);
    fileErrors(noStartTime & ~noStopTime) = {noStartTimeStr};
    fileErrors(noStopTime & ~noStartTime) = {noStopTimeStr};
    fileErrors(noStopTime & noStartTime) = {[noStartTimeStr, '; ', noStopTimeStr]};
    
    % create bad file table
    isBadFile = ~strcmp(fileErrors, '');
    badFiles = table(wavPaths(isBadFile), fileErrors(isBadFile), 'VariableNames',{'Path','Error'});
    
    % compile ranges of valid files
    dtRecs = [dtStarts(~isBadFile), dtStops(~isBadFile)];
    
    % compile start and end times of gaps between files
    dtGaps = [dtRecs(1:(end-1),2), dtRecs(2:end,1)];
    gapDiff = dtGaps(:,2) - dtGaps(:,1);
    
    % issue warning if negative gaps are detected
    if any(gapDiff < seconds(0))
        warning('The recording time series does not appear to be monotonic. Measured gaps may be inaccurate.')
    end
    
    % limit gaps to those within tolerance
    dtGaps = dtGaps(gapDiff >= tol, :);
    
    % plot recording areas
    if doPlot
        
        % create vectors for full timeseries
        xFull = repelem([dtRecs(1,1);dtRecs(end,end)], 2);
        yFull = [0;1;1;0];
        
        % create vectors for periods where there are no gaps
        xNoGaps = repelem([dtRecs(1,1);reshape(dtGaps',numel(dtGaps),1);dtRecs(end,end)], 2);
        yNoGaps = repmat([0;1;1;0], size(dtGaps,1)+1, 1);
        
        % make figure
        f = figure();
        ax = axes('Parent',f);
        ax.NextPlot = 'add';
        area(ax, xFull, yFull, 'FaceColor',[0.6,0.6,0.6], 'EdgeColor','none');
        area(ax, xNoGaps, yNoGaps, 'FaceColor',lines(1), 'EdgeColor','none');
        
        ax.YTick = [];
        xlabel(ax, 'Date-Time')
        title(ax, 'Recording Coverage')
        %zoom(f, 'xon')
    end
end