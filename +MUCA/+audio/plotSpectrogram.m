function varargout = plotSpectrogram(varargin)
%
% Plots a spectrogram from a section of an input audio file.
% If a particular section is not specified, will plot the full file.
%
% SYNTAX:
%   plotSpectrogram(filepath)
%   plotSpectrogram(filepath, Name, Value)
%   plotSpectrogram(ax, __)
%   ax = plotSpectrogram(__)
%
% INPUT ARGUMENTS:
%   Required
%   .......................................................................
%   "filepath" - Char string representing path to an audio file
%   .......................................................................
%
%   Optional (positional arguments)
%   .......................................................................
%   "ax" - Handle of a valid Axes object. If not specified, will use the
%       current Axes (or create a new one if there are no existing Axes).
%   .......................................................................
%
%   Optional (Name-Value pairs)
%   .......................................................................
%   "Channel" - Integer specifying the audio file channel to use, if the 
%       file has multiple channels. Default is 1.
%   .......................................................................
%   "SampleRange" - Interval to extract audio data from, in samples. Must
%       be a 2-element vector of monotonic integers. "SampleRange" cannot
%       be specified at the same time as "TimeRange".
%   .......................................................................
%   "TimeRange" - Time interval to extract audio data from. May be
%       specified as either numeric seconds, a duration object, or (if the
%       audio file contains a timestamp) a datetime object. All cases must
%       be entered as a 2-element vector of increasing values. "TimeRange"
%       cannot be specified at the same time as "SampleRange".
%   .......................................................................
%   "FreqRange" - Range of frequencies to display on the spectrogram. Must
%       be specified as a 2-element vector of increasing numeric values.
%       Default is [0, fs/2].
%   .......................................................................
%   "SpecParams" - Spectrogram parameters. This may be either a char string
%       specifying one of the available pre-determined sets of parameters, 
%       or a struct with custom parameters. Available presets include:
%           'meridian' (default), 'vanderlaan', 'flogeras', 'baumgartner',
%           'pamlab_general', and 'pamlab_whistles'.
%       If using a struct, it must contain exactly the following fields:
%           'win' - window size
%           'ovl' - amount of overlapping samples
%           'nfft' - samples to use in FFT
%   .......................................................................
%   "Smooth" - Boolean value (true/false) specifying whether or not the 
%       spectrogram should be smoothed with a Gaussian kernal. Default is
%       false.
%   .......................................................................
%   "LogFreqs" - Boolean value (true/false) specifying whether or not the 
%       frequency axis should use a logarithmic scale. Default is false.
%   .......................................................................
%   "ColMap" - Char string specifying the colormap to use, or a matrix of 
%       colour values. Supported strings are identical to those accepted by
%       MATLAB's "colormap" function. Default is 'parula'.
%           TIP: To invert the scale of a colormap, input the matrix
%           directly using flipud. For example:
%               flipud(colormap('gray'))
%   .......................................................................
%   "XAxisScale" - Char string specifying the scale of the X axis. The
%       available options are:
%           'abstime' - plot absolute date-times. Only works if the audio
%               file has a timestamp.
%           'reltime' - plot relative time (using durations), where the 
%               audio file starts at 00:00:00.
%           'sec' - plot relative time in numeric seconds, where the audio 
%               file starts at 0.
%           'auto' - automatically choose between 'abstime' and 'reltime'.
%               This will prioritize 'abstime' whenever possible. This is 
%               the default.
%   .......................................................................
%   "XAxisStart" - Overrides the start time of the spectrogram plot (i.e., 
%       the time at x = 0). This only works if using relative time. May be 
%       specified as either a duration object or numeric seconds.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   Optional
%   .......................................................................
%   "ax" - Handle of Axes object on which spectrogram is plotted.
%   .......................................................................
%
% OUTPUT FILES:
%   <none>
%
% DEPENDENCIES:
%   MUCA.time.readDateTime
%   
%
% Written by Wilfried Beslin
% Last updated 2023-11-30 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import MUCA.time.readDateTime

    nargoutchk(0,1)

    % parse 1st input argument, check if it's an axis or not
    arg1 = varargin{1};
    if isa(arg1, 'matlab.graphics.axis.Axes')
        ax = arg1;
        args = varargin(2:end);
    else
        ax = [];
        args = varargin;
    end

    % parse remaining arguments
    p = inputParser;

    p.addRequired('filepath', @isfile)
    p.addParameter('Channel', 1, @(a) validateattributes(a,{'numeric'},{'positive','integer','scalar'}))
    p.addParameter('SampleRange', [], @(a) validateattributes(a,{'numeric'},{'positive','integer','increasing','numel',2}))
    p.addParameter('TimeRange', [], @(a) validateattributes(a,{'numeric','duration','datetime'},{'numel',2}))
    p.addParameter('FreqRange', [], @(a) validateattributes(a,{'numeric'},{'nonnegative','increasing','numel',2}))
    p.addParameter('SpecParams', 'meridian', @(a) validateattributes(a,{'char','struct'},{}))
    p.addParameter('Smooth', false, @(a) validateattributes(a,{'logical'},{'scalar'}))
    p.addParameter('LogFreqs', false, @(a) validateattributes(a,{'logical'},{'scalar'}))
    p.addParameter('ColMap', 'parula', @(a) validateattributes(a,{'char','numeric'},{})) % may be either a colormap string or matrix
    p.addParameter('XAxisScale', 'auto', @(a) validateattributes(a,{'char'},{'row'})) % options are 'auto', 'abstime', or 'reltime'
    p.addParameter('XAxisStart', [], @(a) validateattributes(a,{'numeric','duration'},{'scalar'}))

    p.parse(args{:})

    filepath = p.Results.filepath;
    [~, filename, ~] = fileparts(filepath);
    channel = p.Results.Channel;
    if ~isempty(p.Results.SampleRange)
        assert(isempty(p.Results.TimeRange), 'Cannot specify both ''SampleRange'' and ''TimeRange''')
        interval = p.Results.SampleRange;
        intervalType = 'samples';
    elseif ~isempty(p.Results.TimeRange)
        interval = p.Results.TimeRange;
        intervalType = 'time';
    else
        interval = [];
        intervalType = 'none';
    end
    fRange = p.Results.FreqRange;
    specParamOpt = p.Results.SpecParams;
    doSmooth = p.Results.Smooth;
    doLogFreqs = p.Results.LogFreqs;
    colmap = p.Results.ColMap;
    xAxisScaleOpt = lower(p.Results.XAxisScale);
    xAxisStart = p.Results.XAxisStart; % can be duration or numeric
    
    % get file timestamp if there is one
    fileTimestamp = readDateTime(filename);

    % get sampling rate and number of samples of input file
    filedata = audioinfo(filepath);
    fs = filedata.SampleRate;
    numSamplesTotal = filedata.TotalSamples;

    % process spectrogram parameters
    if ischar(specParamOpt)
        specParams = getPresetSpecParams(specParamOpt, fs);
    else
        specParams = specParamOpt;
    end
    validateSpecParams(specParams)
    sampleStep = specParams.win - specParams.ovl;

    % get sample range to isolate
    sampleBuffer = specParams.win; % spectrogram edges may not show if this is too small
    timeBuffer = duration(0, 0, sampleBuffer./fs);
    switch intervalType
        case 'samples'
            sampleRange = interval;
        case 'time'
            sampleRange = processTimeRange(interval, fs, fileTimestamp);
        case 'none'
            sampleRange = [1, numSamplesTotal];
    end
    extendedSampleRange = sampleRange + [-sampleBuffer, sampleBuffer];
    timeRange = duration(0, 0, (sampleRange-1)./fs);
    timeSpan = diff(timeRange);
    
    % keep track of regions that are out of bounds
    excessLeft = max([0, 1 - extendedSampleRange(1)]);
    excessRight = max([0, extendedSampleRange(end) - numSamplesTotal]);
    inputSampleRange = extendedSampleRange + [excessLeft, -excessRight];
    
    % read audio data
    x = audioread(filepath, inputSampleRange);
    x = x(:, channel);
    
    % pad with zeros if needed
    x = [zeros(excessLeft,1); x; zeros(excessRight,1)];
    
    % compute spectrogram
    [~, specF, specT, specP] = spectrogram(x, specParams.win, specParams.ovl, specParams.nfft, fs);
    %%% Note that P here is the PSD, using the default 'spectrumType' 
    %%% argument
    
    % According to the MATLAB doc:
    %   "The time values in <specT> correspond to the midpoint of each
    %   segment."
    % However, the surf() function does not plot bins - the data points
    % correspond to edges. Therefore, create a new variable that 
    % corresponds to the start of each segment instead. 
    % Also, shift the time vector based on buffer.
    specTAdjusted = specT - (sampleStep./fs)/2;
    specTAdjusted = specTAdjusted - seconds(timeBuffer);
    
    % limit data to selected frequency range
    if isempty(fRange)
        fRange = [0, fs/2];
    end
    fStep = fs./specParams.nfft;
    fInRange = specF >= fRange(1) - fStep & specF <= fRange(end) + fStep;
    specF = specF(fInRange);
    specP = specP(fInRange,:);

    % smoothing
    if doSmooth
        smoothKernel = [1,2,1; 2,4,2; 1,2,1];
        specP = conv2(specP, smoothKernel);
        specP = specP(2:(end-1), 2:(end-1));
    end
    
    % change to dB scale
    specPLog = 10*log10(specP);
    
    % normalize to range [0,1]
    %%% The other normalization method is to use mean of 0 and SD of 1
    %specPLogNorm = specPLog - min(specPLog(:));
    %specPLogNorm = specPLog/max(specPLog(:));
    
    % set x-axis scale for plot
    if strcmp(xAxisScaleOpt, 'auto')
        % check if file has a timestamp; if so, use absolute time,
        % otherwise use relative time.
        if isnat(fileTimestamp)
            xAxisScale = 'reltime';
        else
            xAxisScale = 'abstime';
        end
    else
        xAxisScale = xAxisScaleOpt;
    end
    
    %relTStart = duration(0, 0, (extendedSampleRange(1)-1)./fs);
    relTStart = timeRange(1);
    switch xAxisScale
        case 'abstime'
            if ~isempty(xAxisStart)
                warning('Ignoring ''XAxisStart'' since absolute time will be used')
            end
            tStartPlot = relTStart + fileTimestamp;
            specTPlot = duration(0, 0, specTAdjusted) + tStartPlot;
            tRangePlot = [tStartPlot, tStartPlot + timeSpan];
            xLabelStr = 'Date-Time';
        case 'reltime'
            if ~isempty(xAxisStart)
                if isduration(xAxisStart)
                    tStartPlot = xAxisStart;
                elseif isnumeric(xAxisStart)
                    tStartPlot = duration(0, 0, xAxisStart);
                end
            else
                tStartPlot = relTStart;
            end
            specTPlot = duration(0, 0, specTAdjusted) + tStartPlot;
            tRangePlot = [tStartPlot, tStartPlot + timeSpan];
            xLabelStr = 'Time [hh:mm:ss]';
        case 'sec'
            % display as seconds
            if ~isempty(xAxisStart)
                if isduration(xAxisStart)
                    tStartPlot = seconds(xAxisStart);
                elseif isnumeric(xAxisStart)
                    tStartPlot = xAxisStart;
                end
            else
                tStartPlot = relTStart;
            end
            specTPlot = specTAdjusted + tStartPlot;
            tRangePlot = [tStartPlot, tStartPlot + seconds(timeSpan)];
            xLabelStr = 'Time [s]';
        otherwise
            error('''%s'' is not a supported option for XAxisScale', xAxisScale)
    end

    % initialize plot
    if isempty(ax)
        ax = gca();
    end
    cla(ax);
    ax.NextPlot = 'add';

    % do plot
    surf(ax, specTPlot, specF, specPLog, 'EdgeColor', 'none');
    axis(ax, 'xy');
    view(ax, 0,90);
    colormap(ax, colmap);
    
    % tweak plot
    
    %%% log frequency scale
    if doLogFreqs
        ax.YScale = 'log';
    end
    
    %%% set y axis label and limits
    ylabel(ax, 'Frequency [Hz]');
    ylim(ax, fRange)
    
    %%% set x axis label and limits
    xlabel(ax, xLabelStr)
    xlim(ax, tRangePlot);
    
    %%% add box
    box(ax, 'on')
    
    %%% show tick marks above plot
    set(ax, 'Layer', 'Top')
    
    % set output
    if nargout > 0
        varargout{1} = ax;
    end

end


%% getPresetSpecParams ----------------------------------------------------
function params = getPresetSpecParams(presetName, fs)
% Get pre-defined spectrogram parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch presetName
        case 'meridian'
            specparam_win = round(fs*0.256);
            specparam_ovl = specparam_win - round(fs*0.032);
            specparam_nfft = 2^nextpow2(fs*(2048/8000));
        case 'vanderlaan'
            specparam_win = round(fs*0.2);
            specparam_ovl = specparam_win - round(fs*0.064);
            specparam_nfft = 2^nextpow2(fs/5);
        case 'flogeras'
            specparam_win = 256; % scipy default
            specparam_ovl = round(0.5*specparam_win);
            specparam_nfft = 1024;
        case 'baumgartner'
            specparam_win = round(fs*(640/2048));
            specparam_ovl = round(0.8*specparam_win);
            specparam_nfft = specparam_win;
        case 'pamlab_general'
            specparam_win = round(fs*0.128);
            specparam_ovl = specparam_win - round(fs*0.032);
            specparam_nfft = 2^nextpow2(fs/2);
        case 'pamlab_whistles'
            specparam_win = round(fs*0.05);
            specparam_ovl = specparam_win - round(fs*0.01);
            specparam_nfft = 2^nextpow2(fs/4);
        otherwise
            error('preset spectrogram parameters ''%s'' is not recognized', presetName);
    end
    
    params = struct(...
        'win', specparam_win,...
        'ovl', specparam_ovl,...
        'nfft', specparam_nfft);
end


%% validateSpecParams -----------------------------------------------------
function validateSpecParams(params)
% Make sure spectrogram parameters are in an acceptable format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assert(isstruct(params), '''SpecParams'' must be either a predefined string or a struct with fields ''win'', ''ovl'', and ''nfft''');
    expectedFields = {'win'; 'ovl'; 'nfft'};
    structFields = fieldnames(params);
    assert(numel(structFields) == numel(intersect(expectedFields, structFields)), 'Invalid spectrogram parameter fields. Fields must be ''win'', ''ovl'', and ''nfft''.')
end


%% processTimeRange -------------------------------------------------------
function sampleRange = processTimeRange(interval, fs, dtStart)
% Convert absolute or relative time range into samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % process absolute or relative time forms
    if isdatetime(interval)
        % absolute time; use file timestamp
        if isnat(dtStart)
            error('Cannot process absolute date-time range because timestamp could not be read')
        end
        intSeconds = seconds(interval - dtStart);
    elseif isduration(interval)
        % relative time as duration object
        intSeconds = seconds(interval);
    elseif isnumeric(interval)
        % relative time as numeric seconds
        intSeconds = interval;
    else
        error('Time interval is not in a valid format')
    end
    
    % convert seconds to samples
    sampleRange = intSeconds.*fs + 1;
    
end