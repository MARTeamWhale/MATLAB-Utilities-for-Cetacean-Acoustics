classdef Spectrogram
% Value class to represent an audio spectrogram. Includes time, frequency, 
% and PSD data, as well as specialized plotting methods.
%
% This class tries to keep things simple and does not store spectrogram 
% parameters; it only stores the time vector, frequency vector, PSD matrix,
% and other variables related to or derived from these. That said, the
% class does define a library of preset (sampling frequency-dependent)
% spectrogram parameters that is accessible from the static method
% "get_preset_parameters". Spectrogram objects may be constructed using one
% of the presets, or with a custom set of parameters.
%
% As a convenience feature, this class allows users to edit the PSD matrix
% directly (using both linear and dB scale). This can be useful, e.g., to
% translate the whole matrix such that the maximum equals 0 dB, or to
% normalize frequency bins, or to set a group of values to zero, etc.
%
% Written by Wilfried Beslin
% Last Updated 2025/03/11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Access = public, Dependent)
        psd         % power spectral density matrix (f-by-t)
        psd_db      % power spectral density matrix, in dB (f-by-t)
    end
    properties (SetAccess = private)
        f           % DFT frequency points (column vector), in Hz
    end
    properties (SetAccess = private, Dependent)
        t           % STFT time points (row vector) - may be numeric seconds, duration, or datetime depending on 't_start'. A time point corresponds to the midpoint of a STFT frame.
    end
    properties (SetAccess = private)
        t_relative  % relative time bins (i.e., where first time bin follows 0), in seconds
    end
    properties (Access = public)
        t_start     % start time - may be numeric seconds, duration, or datetime (this will determine the type of 't')
    end
    properties (SetAccess = private)
        df          % frequency step, in Hz
        dt          % time step, in seconds
    end
    properties (SetAccess = private, Dependent)
        nf          % number of frequency bins
        nt          % number of time bins
    end
    properties (Access = private)
        psd_priv    % private instance of power spectral density matrix
    end
    
    
    % Constructor ---------------------------------------------------------
    methods
        function obj = Spectrogram(x, fs, varargin)
        % Creates a Spectrogram object.
        %
        % INPUT PARAMETERS:
        %   Required
        %   ...............................................................
        %   "x" - time series
        %   ...............................................................
        %   "fs" - sampling rate
        %   ...............................................................
        %   
        %   Optional (Name-Value Pairs)
        %   ...............................................................
        %   "params" - may be either a string specifying preset spectrogram
        %       parameters to use, or a struct with the following fields:
        %           'winfcn' - may be either a char string defining a
        %               window function to use (such as 'hamming'), or a
        %               numeric vector that defines a custom window
        %               function (normally with values between 0 and 1).
        %           'winsize' (CONDITIONAL) - an integer that defines the
        %               size of the window (a.k.a. frame), in samples. If
        %               'winfcn' is a vector, this parameter is ignored.
        %           'winparams' (OPTIONAL) - a cell array of arguments for
        %               certain window functions that can accept additional
        %               inputs other than size (for example, the 'blackman'
        %               window accepts a flag for either symmetric or
        %               periodic sampling). This field is completely
        %               optional, and is ignored if 'winfcn' is a vector.
        %           'novl' - The number of samples to overlap between
        %               successive STFT frames. This must be less than
        %               'winsize' (or the length of 'winfcn' if 'winfcn' is
        %               a vector).
        %           'nfft' - The number of points to use in FFT. It must be
        %               greater than or equal to 'winsize' (or the length
        %               of 'winfcn' if 'winfcn' is a vector).
        %       For string inputs, the valid options are defined in the
        %       static method "get_preset_parameters" at the end of this
        %       class file.
        %       The default is 'meridian'.
        %   ...............................................................
        %   "smooth" - Boolean value specifying whether or not the
        %       spectrogram should be smoothed with a Gaussian kernal.
        %       Default is false.
        %   ...............................................................
        %   "t_start" - Reference start time for t vector. May be either
        %       numeric seconds, a duration object, or a datetime object.
        %       The data type of t_start determines the data type of the
        %       whole t vector.
        %
        %       NOTE: since t corresponds to the midpoint of STFT frames,
        %       t(1) is not equal to t_start. Rather, all values of t are
        %       measured relative to t_start.
        %
        %       Default is 0 (numeric).
        %   ...............................................................
        %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            p = inputParser;
            p.addRequired('x', @(v)validateattributes(v,{'numeric'},{'vector'}))
            p.addRequired('fs', @(v)validateattributes(v,{'numeric'},{'positive','integer','scalar'}))
            p.addParameter('params', 'meridian');
            p.addParameter('smooth', false);
            p.addParameter('t_start', 0);

            p.parse(x, fs, varargin{:});
            params_input = p.Results.params;
            do_smooth = p.Results.smooth;
            t0 = p.Results.t_start;

            % process spectrogram parameters
            if isstruct(params_input)
                params_struct = params_input;
            else
                params_struct = obj.get_preset_parameters(params_input, fs);
            end
            [winvec, novl, nfft] = obj.parse_parameters(params_struct);

            % compute spectrogram
            %%% using the default 'spectrumType' argument, which returns
            %%% the PSD
            [~, spec_f, spec_t, spec_psd] = spectrogram(x, winvec, novl, nfft, fs);
            
            % determine time and frequency step
            delta_t = (numel(winvec) - novl)/fs;
            delta_f = fs/nfft;
            
            % set object properties
            obj.psd_priv = spec_psd;
            obj.f = spec_f;
            obj.t_relative = spec_t;
            obj.t_start = t0;
            obj.df = delta_f;
            obj.dt = delta_t;

            % smoothing
            if do_smooth
                obj = obj.smooth();
            end
        end
    end
    
    
    %% "SET" METHODS ======================================================
    methods
        % set.psd .........................................................
        function obj = set.psd(obj, val)
            assert(all(size(val) == [obj.nf, obj.nt]))
            obj.psd_priv = val;
        end
        
        % set.psd_db ......................................................
        function obj = set.psd_db(obj, val)
            assert(all(size(val) == [obj.nf, obj.nt]))
            obj.psd_priv = 10.^(val./10);
        end
        
        % set.t_start .....................................................
        function obj = set.t_start(obj, val)
            validateattributes(val, {'numeric','duration','datetime'}, {'scalar'})
            obj.t_start = val;
        end
    end


    %% "GET" METHODS ======================================================
    methods
        % get.psd .........................................................
        function val = get.psd(obj)
            val = obj.psd_priv;
        end
        
        % get.psd_db ......................................................
        function val = get.psd_db(obj)
            val = 10*log10(obj.psd_priv);
        end
        
        % get.t ...........................................................
        function val = get.t(obj)
            t0 = obj.t_start;
            if isnumeric(t0)
                t_base = obj.t_relative;
            else
                t_base = seconds(obj.t_relative);
            end
            val = t0 + t_base; % note: for duration objects, the first variable in this operation determines the format of the output.
        end
        
        % get.nf ..........................................................
        function val = get.nf(obj)
            val = numel(obj.f);
        end
        
        % get.nt ..........................................................
        function val = get.nt(obj)
            val = numel(obj.t_relative);
        end
    end
    
    
    %% PUBLIC METHODS =====================================================
    methods (Access = public)
        % smooth ----------------------------------------------------------
        function obj = smooth(obj)
        % Smooth spectrogram with a 2D Gaussian kernel function.
        % The method is based on:
        %   Baumgartner and Mussoline (2011), "A generalized baleen whale
        %   call detection and classification system". J. Acoust. Soc. Am. 
        %   129: 2889-2902. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % define Gaussian smoothing kernel
            M = [1,2,1; 2,4,2; 1,2,1];

            % duplicate the time and frequency edges of the spectrogram to
            % prepare for edge effects
            psd_anal = obj.psd;
            psd_anal = [psd_anal(:,1), psd_anal, psd_anal(:,end)];
            psd_anal = [psd_anal(1,:); psd_anal; psd_anal(end,:)];

            % convolve the spectrogram matrix with the smoothing kernel
            %%% the matrix must also be divided by the sum of the smoothing
            %%% kernel elements to preserve the scale
            psd_smooth = conv2(psd_anal, M)./sum(M(:));

            % remove the edges introduced by both edge duplication and
            % convolution (2 levels)
            psd_smooth = psd_smooth(3:(end-2), 3:(end-2));

            % set smoothed matrix as the PSD of the new object
            obj.psd = psd_smooth;
        end
        
        % truncate --------------------------------------------------------
        function obj = truncate(obj, varargin)
        % Truncate a spectrogram in time and/or frequency.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            p = inputParser;
            p.addParameter('t_min', []);
            p.addParameter('t_max', []);
            p.addParameter('f_min', []);
            p.addParameter('f_max', []);
            
            p.parse(varargin{:});
            
            if ~all(ismember({'t_min','t_max','f_min','f_max'}, p.UsingDefaults))
                % process t_min
                t_min_old = obj.t(1);
                t_min_new = p.Results.t_min;
                if isempty(t_min_new)
                    t_min_new = t_min_old;
                elseif t_min_new < t_min_old
                    warning('Ignoring out-of-bounds value: t_min')
                    t_min_new = t_min_old;
                end

                % process t_max
                t_max_old = obj.t(end);
                t_max_new = p.Results.t_max;
                if isempty(t_max_new)
                    t_max_new = t_max_old;
                elseif t_max_new > t_max_old
                    warning('Ignoring out-of-bounds value: t_max')
                    t_max_new = t_max_old;
                end

                % make sure times are within range
                assert(t_min_new < t_max_new, 'Non-monotonic time range')

                % process f_min
                f_min_old = obj.f(1);
                f_min_new = p.Results.f_min;
                if isempty(f_min_new)
                    f_min_new = f_min_old;
                elseif f_min_new < f_min_old
                    warning('Ignoring out-of-bounds value: f_min')
                    f_min_new = f_min_old;
                end

                % process f_max
                f_max_old = obj.f(end);
                f_max_new = p.Results.f_max;
                if isempty(f_max_new)
                    f_max_new = f_max_old;
                elseif f_max_new > f_max_old
                    warning('Ignoring out-of-bounds value: f_max')
                    f_max_new = f_max_old;
                end

                % make sure frequencies are within range
                assert(f_min_new < f_max_new, 'Non-monotonic frequency range')

                % calculate truncated ranges
                t_keep = obj.t >= t_min_new & obj.t <= t_max_new;
                f_keep = obj.f >= f_min_new & obj.f <= f_max_new;
                
                % determine new relative time vector
                t_rel_new = obj.t_relative(1:sum(t_keep));
                
                % determine start time shift if minimum time was changed
                %%% But only if t_start is something other than 0
                if ~t_keep(1) && (isdatetime(obj.t_start) || obj.t_start ~= 0)
                    t1_old = obj.t(1);
                    t1_new = obj.t(find(t_keep, 1, 'first'));
                    delta_t1 = t1_new - t1_old;
                    t_start_new = obj.t_start + delta_t1;
                else
                    t_start_new = obj.t_start;
                end

                % apply truncation
                obj.t_start = t_start_new;
                obj.t_relative = t_rel_new;
                obj.f = obj.f(f_keep);
                obj.psd = obj.psd(f_keep, t_keep);
            end
        end
        
        % plot ------------------------------------------------------------
        function varargout = plot(obj, varargin)
        % Produce a plot of the spectrogram.
        %
        % SYNTAX:
        %   plot(obj)
        %   plot(obj, ax)
        %   plot(__, Name,Value)
        %   ax = plot(__)
        %
        % INPUT PARAMETERS:
        %   Optional (positional arguments)
        %   ...............................................................
        %   "ax" - Axes object. If none, uses current Axes or creates one.
        %   ...............................................................
        %   
        %   Optional (Name-Value Pairs)
        %   ...............................................................
        %   "log_freqs" - Boolean value, determines if frequencies are
        %       plotted on log scale or not
        %   ...............................................................
        %   The rest are identical to the Name-Value pairs supported by the
        %   "surf" method in this class. This includes the "type"
        %   parameter, which may be either 'pixels' (default) or 'smooth';
        %   the remaining args are passed to the "truncate" method.
        %   ...............................................................
        %
        % OUTPUT ARGUMENTS:
        %   ...............................................................
        %   "ax" - Axes object
        %   ...............................................................
        %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nargoutchk(0, 1)
            
            % parse input
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('ax', gca)
            p.addParameter('log_freqs', false)
                
            p.parse(varargin{:});
            ax = p.Results.ax;
            log_freqs = p.Results.log_freqs;
            surf_args = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
        
            % make plot
            obj.surf(ax, surf_args{:});
            
            % tweak plot
            axis(ax, 'xy');
            view(ax, 0, 90);
            if log_freqs
                ax.YScale = 'log';
            end
            xlabel(ax, 'Time [s]')
            ylabel(ax, 'Frequency [Hz]');
            box(ax, 'on');
            grid(ax, 'off');
            set(ax, 'Layer', 'Top'); % show tick marks above plot
            
            % set axes limits
            axis(ax, 'tight');
            
            % return plot if needed
            if nargout > 0
                varargout = {ax};
            end
        end
        
        % surf ------------------------------------------------------------
        function varargout = surf(obj, varargin)
        % create a Surface object from spectrogram data
        %
        % SYNTAX:
        %   surf(obj)
        %   surf(obj, ax)
        %   surf(__, Name,Value)
        %   s = surf(__)
        %
        % INPUT PARAMETERS:
        %   Optional (positional arguments)
        %   ...............................................................
        %   "ax" - Axes object. If none, uses current Axes or creates one.
        %   ...............................................................
        %   
        %   Optional (Name-Value Pairs)
        %   ...............................................................
        %   "type" - char string that describes the type of surface plot to
        %       make to represent the spectrogram. The options are:
        %
        %       'pixels' (default) - each face of the surface is centred on
        %           a specific time-frequency point, and has solid
        %           colouring that corresponds to the PSD value of that
        %           point.
        %       'smooth' - faces represent transitions between
        %           time-frequency points, as opposed to the points
        %           themselves. The colour of each face is non-uniform, and
        %           changes gradually from one vertex to another.
        %   ...............................................................
        %   The rest are identical to the "truncate" method in this class
        %   ...............................................................
        %
        % OUTPUT ARGUMENTS:
        %   ...............................................................
        %   "s" - Surface object
        %   ...............................................................
        %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            nargoutchk(0, 1)

            % parse input
            p = inputParser;
            p.KeepUnmatched = true;
            p.addOptional('ax', gca)
            p.addParameter('type', 'pixels')

            p.parse(varargin{:});
            ax = p.Results.ax;
            surf_type = lower(p.Results.type);
            surf_type = validatestring(surf_type, {'pixels','smooth'});
            trunc_args = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
            
            % truncate spectrogram if needed
            obj_plot = obj.truncate(trunc_args{:});
            
            % prepare data for plotting, depending on surface type
            switch surf_type
                case 'pixels'
                    % For this surface type, we want faces that are centred
                    % on time-frequency points. Remember that the data
                    % points in a Surface object correspond to vertices,
                    % not faces; thus, to centre each face, we need to
                    % shift the time points by -dt/2. Likewise, the
                    % frequency points must also be shifted by -df/2.
                    % Finally, we also need to add an additional time point
                    % and an additional frequency point at the end of each
                    % vector to be able to generate faces for the last time
                    % and frequency bins, respectively.
                    %
                    % This type of surface will have a 'FaceColor' value
                    % set to 'flat', which means that each face will be
                    % coloured according to the value of its lower-left
                    % vertex (i.e., the time-frequency point that the face
                    % represents).
                    face_col_type = 'flat';

                    if isnumeric(obj.t)
                        delta_t = obj_plot.dt;
                    else
                        delta_t = seconds(obj_plot.dt);
                    end
                    
                    %%% shift time points and add extra one at the end
                    t_plot = obj_plot.t - delta_t/2;
                    t_plot = [t_plot, t_plot(end) + delta_t];
                    
                    %%% shift frequency points and add extra one at the end
                    f_plot = obj_plot.f - obj_plot.df/2;
                    f_plot = [f_plot; f_plot(end) + obj_plot.df];
                    
                    %%% add an additional row and column at the ends of the
                    %%% PSD matrix
                    psd_plot = obj_plot.psd_db;
                    psd_plot = [psd_plot, psd_plot(:,end)];
                    psd_plot = [psd_plot; psd_plot(end,:)];
                case 'smooth'
                    % For this surface type, the surface data points will
                    % correspond exactly to the spectrogram time-frequency
                    % points. The faces will represent transitions between
                    % these points; but in order to make the spectrogram
                    % colours properly correspond to the power at each
                    % point, we must set the 'FaceColor' property to
                    % 'interp'.
                    face_col_type = 'interp';

                    t_plot = obj_plot.t;
                    f_plot = obj_plot.f;
                    psd_plot = obj_plot.psd_db;
                otherwise
                    error('Invalid surface type')
            end
            
            % create surface
            s = surf(ax, t_plot, f_plot, psd_plot, 'EdgeColor','none', 'FaceColor',face_col_type);
            
            % return surface object if needed
            if nargout > 0
                varargout = {s};
            end
        end
    end


    %% STATIC METHODS (Public) ============================================
    methods (Access = public, Static)
        % parse_parameters ------------------------------------------------
        function [winvec, novl, nfft] = parse_parameters(params_struct)
        % Takes a parameter struct as input and returns variables to be
        % used as input for the "spectrogram" function. The rules of the
        % parameter struct are explained in the header of the class
        % constructor.
        % This method is really intended for internal use only.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % make sure the fields are all expected
            valid_param_fields = {...
                'winfcn',...
                'winsize',...
                'winparams',...
                'novl',...
                'nfft'}; 
            param_fields = fieldnames(params_struct);
            assert(all(ismember(param_fields, valid_param_fields)), 'The parameter struct contains invalid fields')

            % check if 'winfcn' is a string or a vector, and process
            % accordingly
            if ischar(params_struct.winfcn) || isstring(params_struct.winfcn)
                % make sure 'winsize' is specified
                assert(ismember('winsize',param_fields), 'The field ''winsize'' must be specified when ''winfcn'' is a string');

                % process window string input
                winfcn_name = lower(char(params_struct.winfcn));

                %%% convert name to a MATLAB function handle
                winfcn_fcn = str2func(winfcn_name);

                %%% extract any optional window arguments
                if ismember('winparams', param_fields)
                    assert(iscell(params_struct.winparams), 'Expected ''winparams'' to be a cell array')
                    winfcn_args = params_struct.winparams;
                else
                    winfcn_args = {};
                end

                %%% create the window vector
                winvec = window(winfcn_fcn, params_struct.winsize, winfcn_args{:});
            else
                % process window vector input
                %%% make sure the window is a numeric vector
                validateattributes(params_struct.winfcn, {'numeric'},{'vector'});

                %%% if 'winsize' or 'winparams' were specified, issue a
                %%% warning that they will be ignored
                if any(ismember({'winsize','winparams'}, param_fields))
                    warning('''winsize'' and/or ''winparams'' exist, but ''winvec'' was specified directly as a vector and does not require these parameters. They will be ignored.')
                end
                
                %%% assign the window vector
                winvec = params_struct.winfcn;
            end

            % validate 'novl'
            validateattributes(params_struct.novl, {'numeric'}, {'integer','scalar','nonnegative','<',numel(winvec)}, '', 'novl');
            novl = params_struct.novl;

            % validate 'nfft'
            validateattributes(params_struct.nfft, {'numeric'}, {'integer','scalar','positive','>=',numel(winvec)}, '', 'nfft');
            nfft = params_struct.nfft;
        end

        % get_preset_parameters -------------------------------------------
        function params = get_preset_parameters(preset_name, fs)
        % Returns a struct with pre-determined spectrogram parameters.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % parse input
            p = inputParser();
            p.addRequired('preset_name', @(v)validateattributes(v,{'char','string'},{}))
            p.addRequired('fs', @(v)validateattributes(v,{'numeric'},{'scalar','positive','integer'}))
            p.parse(preset_name, fs)

            preset_name = lower(char(preset_name));
            valid_presets = {...
                'meridian',...
                'bio_narw_analysis',...
                'bio_muppet',...
                'vanderlaan_2003',...
                'baumgartner_2011',...
                'pamlab_general',...
                'pamlab_longcalls',...
                'pamlab_clicks',...
                'pamlab_hiresclicks',...
                'pamlab_whistles'};
            preset_name = validatestring(preset_name, valid_presets);

            % initialize output
            params = struct(...
                'winfcn', [],...
                'winsize', [],...
                'novl', [],...
                'nfft', []);

            % process parameters
            switch preset_name
                case 'meridian'
                    % parameters used by the MERIDIAN group to process
                    % right whale upcalls, such as in Kirsebom et al.
                    % (2020).
                    params.winfcn = 'hamming'; % not actually sure what it is
                    params.winsize = round(fs*0.256);
                    params.novl = params.winsize - round(fs*0.032);
                    params.nfft = 2^nextpow2(fs*(2048/8000));
                case 'bio_narw_analysis'
                    % parameters typically used at BIO to detect right
                    % whale upcalls with LFDCS
                    params.winfcn = 'hann';
                    params.winsize = round(fs*(512/2000));
                    params.novl = round(0.75*params.winsize);
                    params.nfft = params.winsize;
                case 'bio_muppet'
                    % default parameters used in MUPPET for blue whale call
                    % analysis
                    params.winfcn = 'hamming';
                    params.winsize = round(fs*(8192/8000));
                    params.novl = round(params.winsize*0.9);
                    params.nfft = params.winfcn;
                case 'flogeras_aural'
                    % parameters used by Dave Flogeras's "aural-features"
                    % python package (MIGHT NEED REVISION). 
                    params.winfcn = 'hamming'; % not actually sure what it is
                    params.winsize = 256; % scipy default
                    params.novl = round(0.5*params.winsize);
                    params.nfft = 1024;
                case 'vanderlaan_2003'
                    % parameters used in Vanderlaan et al. (2003) to
                    % characterize right whale calls in the Bay of Fundy.
                    params.winfcn = 'hamming'; % not actually sure what it is
                    params.winsize = round(fs*0.2);
                    params.novl = params.winsize - round(fs*0.064);
                    params.nfft = 2^nextpow2(fs/5);
                case 'baumgartner_2011'
                    % parameters used by Baumgartner and Mussoline (2011)
                    % in their data analysis during their presentation of
                    % LFDCS
                    params.winfcn = 'hann';
                    params.winsize = round(fs*(640/2048));
                    params.novl = round(0.8*params.winsize);
                    params.nfft = params.winsize;
                case 'pamlab_general'
                    % parameters from the "General (all call types)"
                    % single-part setting in PAMlab
                    params.winfcn = 'hamming';
                    params.winsize = round(fs*0.128);
                    params.novl = params.winsize - round(fs*0.032);
                    params.nfft = 2^nextpow2(fs/2);
                case 'pamlab_longcalls'
                    % parameters from the "Long calls (blue whale)"
                    % single-part setting in PAMlab
                    params.winfcn = 'hamming';
                    params.winsize = round(fs*2);
                    params.novl = params.winsize - round(fs*0.5);
                    params.nfft = 2^nextpow2(fs/0.4);
                case 'pamlab_clicks'
                    % parameters from the "Clicks" single-part setting in
                    % PAMlab
                    params.winfcn = 'hamming';
                    params.winsize = round(fs*0.001);
                    params.novl = params.winsize - round(fs*0.0005);
                    params.nfft = 2^nextpow2(fs/128);
                case 'pamlab_hiresclicks'
                    % parameters from the "High-res clicks (beaked whale)"
                    % single-part setting in PAMlab
                    params.winfcn = 'hamming';
                    params.winsize = round(fs*0.000266);
                    params.novl = params.winsize - round(fs*0.00002);
                    params.nfft = 2^nextpow2(fs/512);
                case 'pamlab_whistles'
                    % parameters from the "Whistles" single-part setting in
                    % PAMlab
                    params.winfcn = 'hamming';
                    params.winsize = round(fs*0.05);
                    params.novl = params.winsize - round(fs*0.01);
                    params.nfft = 2^nextpow2(fs/4);
                otherwise
                    error('This preset lacks a definition, please update the code to add one!')
            end
        end
    end
end