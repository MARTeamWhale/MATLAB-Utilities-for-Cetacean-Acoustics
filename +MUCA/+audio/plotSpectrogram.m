function varargout = plotSpectrogram(varargin)
%
% Plots a spectrogram from a section of an input audio file.
% If a particular section is not specified, will plot the full file.
%
% Note that this function relies on the MUCA.dsp.Spectrogram class. Some of
% the input arguments that can be specified here are passed directly to
% MUCA.dsp.Spectrogram.
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
%       audio file contains a timestamp or a starting datetime has been 
%       specified,) a datetime object. All cases must
%       be entered as a 2-element vector of increasing values. "TimeRange"
%       cannot be specified at the same time as "SampleRange".
%   .......................................................................
%   "FreqRange" - Range of frequencies to display on the spectrogram. Must
%       be specified as a 2-element vector of increasing numeric values.
%       Default is [0, fs/2].
%   .......................................................................
%   "SpecParams" - Spectrogram parameters. This may be either a char string
%       specifying one of the available pre-determined sets of parameters
%       in the MUCA.dsp.Spectrogram class, or a struct of parameters
%       accepted by MUCA.dsp.Spectrogram. See the class documentation for
%       details.
%   .......................................................................
%   "Smooth" - Boolean value (true/false) specifying whether or not the 
%       spectrogram should be smoothed with a Gaussian kernal (see
%       MUCA.dsp.Spectrogram for details). Default is false.
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
%               file has a timestamp, or if 'XAxisStart' has been specified
%               as a datetime object.
%           'reltime' - plot relative time (using durations), where the 
%               audio file starts at 00:00:00.
%           'sec' - plot relative time in numeric seconds, where the audio 
%               file starts at 0.
%           'auto' - automatically choose between 'abstime' and 'reltime'.
%               This will prioritize 'abstime' whenever possible. This is 
%               the default.
%   .......................................................................
%   "XAxisStart" - Overrides the start time of the spectrogram plot (i.e., 
%       the time at x = 0). May be specified as either a duration object, 
%       numeric seconds, or a datetime object. Datetime is only supported
%       if 'XAxisScale' is 'abstime' or 'auto'. Specifying 'XAxisStart' 
%       with a datetime object will also override any timestamp that may 
%       be coded in the filename.
%   .......................................................................
%   "CentreFirstSample" - Boolean value (true/false) specifying if the
%       first sample in the range of interest should occur at the centre of
%       a STFT frame (true) or the beginning (false). Setting it at the
%       beginning is the usual behaviour (at least when using MATLAB's
%       "spectrogram" function), but this means that the first few samples
%       may be attenuated by the window function. The default is false.
%   .......................................................................
%   "PlotType" - char string that determines if the spectrogram plot should
%       appear granualar ('pixels') or smooth ('smooth'). This is
%       equivalent to the 'type' parameter in MUCA.dsp.Spectrogram.plot.
%       The default is 'pixels'.
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
%   MUCA.dsp.Spectrogram
%   
%
% Written by Wilfried Beslin
% Last updated 2025-03-12 using MATLAB R2024a
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import MUCA.time.readDateTime
    import MUCA.dsp.Spectrogram

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
    p.addParameter('SampleRange', [], @(a) validateattributes(a,{'numeric'},{'positive','integer','increasing','row','numel',2}))
    p.addParameter('TimeRange', [], @(a) validateattributes(a,{'numeric','duration','datetime'},{'row','numel',2}))
    p.addParameter('FreqRange', [], @(a) validateattributes(a,{'numeric'},{'nonnegative','increasing','row','numel',2}))
    p.addParameter('SpecParams', 'meridian', @(a) validateattributes(a,{'char','struct'},{}))
    p.addParameter('Smooth', false, @(a) validateattributes(a,{'logical'},{'scalar'}))
    p.addParameter('LogFreqs', false, @(a) validateattributes(a,{'logical'},{'scalar'}))
    p.addParameter('ColMap', 'parula', @(a) validateattributes(a,{'char','numeric'},{})) % may be either a colormap string or matrix
    p.addParameter('XAxisScale', 'auto', @(a) validateattributes(a,{'char'},{'row'})) % options are 'auto', 'abstime', or 'reltime'
    p.addParameter('XAxisStart', [], @(a) validateattributes(a,{'numeric','duration','datetime'},{'scalar'}))
    p.addParameter('CentreFirstSample', false, @(a) validateattributes(a,{'logical'},{'scalar'}))
    p.addParameter('PlotType', 'pixels', @(a) validateattributes(a,{'char'},{'row'}))

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
    xAxisStart = p.Results.XAxisStart; % can be duration, numeric, or (conditionally) datetime
    centreFirstSample = p.Results.CentreFirstSample;
    plotType = p.Results.PlotType;
    
    % initialize timestamp (will be obtained later as needed)
    fileTimeStamp = datetime.empty(0,0);
    
    % determine or validate the X-axis scale type
    switch xAxisScaleOpt
        case 'auto'
            % Check if XAxisStart is datetime - if it is, use 'abstime'
            if isdatetime(xAxisStart)
                xAxisScale = 'abstime';
            elseif ~isempty(xAxisStart)
                % if XAxisStart was specified and is something other than
                % datetime, then use 'reltime'
                xAxisScale = 'reltime';
            else
                % check if the file has a timestamp
                fileTimeStamp = readDateTime(filename, 'SuppressWarnings',true);
                if isnat(fileTimeStamp)
                    % no timestamp - use 'reltime'
                    xAxisScale = 'reltime';
                else
                    % timestamp found - use 'abstime'
                    xAxisScale = 'abstime';
                end
            end
            
        case 'abstime'
            % check if XAxisStart was specified - if it was, then it must
            % be datetime
            if ~isempty(xAxisStart)
                assert(isdatetime(xAxisStart), 'XAxisStart must be a "datetime" type if XAxisScale == ''abstime''!')
            else
                % if XAxisStart was not specified, then the file must have
                % a timestamp
                fileTimeStamp = readDateTime(filename, 'SuppressWarnings',true);
                assert(~isnat(fileTimeStamp), sprintf('Could not read timestamp from file "%s".\nA timestamp must exist if XAxisScale == ''abstime'' and XAxisStart is not specified.',filename));
            end
            xAxisScale = xAxisScaleOpt;
            
        case {'reltime','sec'}
            % make sure XAxisStart is not datetime
            assert(~isdatetime(xAxisStart), 'XAxisStart cannot be datetime when XAxisScale == ''reltime'' or ''sec''!')
            xAxisScale = xAxisScaleOpt;
                
        otherwise
            error('''%s'' is not a supported option for XAxisScale', xAxisScaleOpt)
    end

    % get sampling rate and number of samples of input file
    filedata = audioinfo(filepath);
    fs = filedata.SampleRate;
    numSamplesTotal = filedata.TotalSamples;

    % process spectrogram parameters
    %%% while the dsp.Spectrogram constructor could take care of this, the
    %%% parameters must be pre-processed here to ensure that the
    %%% appropriate number of buffer samples required for plotting can be
    %%% determined.
    if ischar(specParamOpt) || isstring(specParamOpt)
        specParams = Spectrogram.get_preset_parameters(specParamOpt, fs);
    else
        specParams = specParamOpt;
    end

    if ischar(specParams.winfcn) || isstring(specParams.winfcn)
        numFrameSamples = specParams.winsize;
    else
        numFrameSamples = numel(specParams.winfcn);
    end
    sampleStep = numFrameSamples - specParams.novl;

    % get the base sample range (i.e. the samples of interest)
    switch intervalType
        case 'samples'
            targetSampleRange = interval;
        case 'time'
            if isempty(fileTimeStamp) && isdatetime(interval)
                fileTimeStamp = readDateTime(filename, 'SuppressWarnings',true);
            end
            targetSampleRange = processTimeRange(interval, fs, fileTimeStamp, xAxisStart);
        case 'none'
            targetSampleRange = [1, numSamplesTotal];
    end
    targetTimeRange = duration(0, 0, (targetSampleRange-1)./fs);
    targetTimeSpan = diff(targetTimeRange);
    numTargetSamples = targetSampleRange(2) - targetSampleRange(1) + 1;

    % determine the number of left-side buffer samples
    %%% this is necessary to ensure that the plotted range does not start
    %%% with a blank. To achieve this, the midpoint of the first frame must
    %%% be equal to or slightly earlier than the index of the first sample.
    %%% But this must be done carefully to make sure that the first sample
    %%% is positioned where we expect it to be in the first full data frame
    %%% (i.e., centred or at the beginning).
    if centreFirstSample
        % in this case, the first STFT frame should start when half (or
        % just over half) of the frame covers the signal
        leftSampleBuffer = floor(numFrameSamples/2);
    else
        % in this case, want the first sample to occur at the beginning of
        % a STFT frame. But there needs to be other frames before this,
        % otherwise the left edge of the plot will be blank. This means
        % there should be a left-side buffer that is equal to enough step
        % sizes that the midpoint of the first frame is equal to or just
        % smaller than 1.
        leftSampleBuffer = sampleStep*ceil((numFrameSamples/2)/sampleStep);
    end
    leftTimeBuffer = duration(0, 0, leftSampleBuffer/fs);

    % determine the number of right-side buffer samples
    %%% similarly to the left side, the right side should contain just
    %%% enough samples to ensure that the midpoint time of the last STFT
    %%% frame is equal to or slightly over the time of the last sample, to
    %%% avoid having a blank on the right side.
    stftInitialFrameStarts = (1 - leftSampleBuffer):sampleStep:numTargetSamples;
    stftInitialFrameMidpoints = stftInitialFrameStarts + numFrameSamples/2;
    numFinalFrames = find(stftInitialFrameMidpoints >= numTargetSamples, 1, 'first');

    stftFrameStarts = stftInitialFrameStarts(1:numFinalFrames);
    stftFrameStops = stftFrameStarts + numFrameSamples;

    rightSampleBuffer = stftFrameStops(end) - numTargetSamples;

    % isolate the part of the time series to use in the STFT
    extendedSampleRange = targetSampleRange + [-leftSampleBuffer, rightSampleBuffer];
    
    % keep track of regions that are out of bounds
    excessLeft = max([0, 1 - extendedSampleRange(1)]);
    excessRight = max([0, extendedSampleRange(end) - numSamplesTotal]);
    inputSampleRange = extendedSampleRange + [excessLeft, -excessRight];
    
    % read audio data
    x = audioread(filepath, inputSampleRange);
    x = x(:, channel);
    
    % pad with zeros if needed
    %** THOUGHT: there could be different padding types: zeros, endpoint
    %** replication, mirroring, random noise... something to experiment
    %** with perhaps.
    x = [zeros(excessLeft,1); x; zeros(excessRight,1)];

    % set time properties based on X-Axis input options
    relTStart = targetTimeRange(1) - leftTimeBuffer;
    switch xAxisScale
        case 'abstime'
            if isempty(xAxisStart)
                % use timestamp
                tSpecStart = relTStart + fileTimeStamp;
            else
                % use specified start time
                tSpecStart = xAxisStart;
            end
            tRangePlot = [tSpecStart, tSpecStart + targetTimeSpan] + leftTimeBuffer;
            xLabelStr = 'Date-Time';
            
        case 'reltime'
            if ~isempty(xAxisStart)
                if isduration(xAxisStart)
                    tSpecStart = xAxisStart;
                elseif isnumeric(xAxisStart)
                    tSpecStart = duration(0, 0, xAxisStart);
                end
            else
                tSpecStart = relTStart;
            end
            tRangePlot = [tSpecStart, tSpecStart + targetTimeSpan] + leftTimeBuffer;
            xLabelStr = 'Time [hh:mm:ss]';
            
        case 'sec'
            % display as numeric seconds
            if ~isempty(xAxisStart)
                if isduration(xAxisStart)
                    tSpecStart = seconds(xAxisStart);
                elseif isnumeric(xAxisStart)
                    tSpecStart = xAxisStart;
                end
            else
                tSpecStart = seconds(relTStart);
            end
            tRangePlot = [tSpecStart, tSpecStart + seconds(targetTimeSpan)] + seconds(leftTimeBuffer);
            xLabelStr = 'Time [s]';
            
        otherwise
            error('''%s'' is not a supported option for XAxisScale', xAxisScale)
    end

    % compute spectrogram
    spec = Spectrogram(x, fs, 'params',specParams, 'smooth',doSmooth, 't_start',tSpecStart);

    % initialize plot
    if isempty(ax)
        ax = gca();
    end
    cla(ax);
    ax.NextPlot = 'add';

    % do plot
    if isempty(fRange)
        fRange = [0, fs/2];
        truncArgs = {};
    else
        truncArgs = {'f_min',fRange(1)-spec.df, 'f_max',fRange(2)+spec.df};
    end
    spec.plot(ax, 'log_freqs',doLogFreqs, 'type',plotType, truncArgs{:})
    colormap(ax, colmap);
    
    % tweak plot
    
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


%% processTimeRange -------------------------------------------------------
function sampleRange = processTimeRange(interval, fs, fileDTStart, xAxisStart)
% Convert absolute or relative time range into samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % process absolute or relative time forms
    if isdatetime(interval)
        % absolute time; use file timestamp
        if isnat(fileDTStart)
            error('Cannot process absolute date-time range because no valid timestamp was found')
        elseif isdatetime(xAxisStart)
            warning('Using absolute time intervals to isolate samples, but a starting date-time for the plot has been specified manually. The times displayed on the plot could be incorrect.')
        end
        intSeconds = seconds(interval - fileDTStart);
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