classdef Spectrogram
% Value class to represent an audio spectrogram. Includes time, frequency, 
% and PSD data, as well as specialized plotting methods.
%
% This class tries to keep things simple and does not store spectrogram 
% parameters; it only stores the time vector, frequency vector, PSD matrix,
% and other variables related to or derived from these.
%
% As a convenience feature, this class allows users to edit the PSD matrix
% directly (using both linear and dB scale). This can be useful, e.g., to
% translate the whole matrix such that the maximum equals 0 dB, or to set a
% group of values to zero, etc.
%
% Written by Wilfried Beslin
% Last Updated 2023/08/19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties (Access = public, Dependent)
        psd         % power spectral density matrix (f-by-t)
        psd_db      % power spectral density matrix, in dB (f-by-t)
    end
    properties (SetAccess = private)
        f           % frequency bins (column vector), in Hz
    end
    properties (SetAccess = private, Dependent)
        t           % time bins (row vector) - may be numeric seconds, duration, or datetime depending on 't_start'
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
    
    
    % Constructor
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
        %       parameters to use, or a struct with fields 'wintype', 
        %       'winsize', 'ovl', and 'nfft'. If using a string, options 
        %       are:
        %           'meridian' (default)
        %           'vanderlaan'
        %           'flogeras'
        %           'baumgartner'.
        %   ...............................................................
        %   "smooth" - Boolean value specifying whether or not the
        %       spectrogram should be smoothed with a Gaussian kernal.
        %       Default is false.
        %   ...............................................................
        %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            p = inputParser;
            p.addParameter('params', 'meridian');
            p.addParameter('smooth', false);
            p.addParameter('t_start', 0); % start time for t vector; determines vector type (numeric, duration, or datetime)

            p.parse(varargin{:});
            params = p.Results.params;
            do_smooth = p.Results.smooth;
            t0 = p.Results.t_start;

            % parameters
            if ischar(params)
                switch params
                    case 'meridian'
                        specparam_winfcn = @hamming; % not actually sure what it is
                        specparam_winsize = round(fs*0.256);
                        specparam_ovl = specparam_winsize - round(fs*0.032);
                        specparam_nfft = 2^nextpow2(fs*(2048/8000));
                    case 'vanderlaan'
                        specparam_winfcn = @hamming; % not actually sure what it is
                        specparam_winsize = round(fs*0.2);
                        specparam_ovl = specparam_winsize - round(fs*0.064);
                        specparam_nfft = 2^nextpow2(fs/5);
                    case 'flogeras'
                        specparam_winfcn = @hamming; % not actually sure what it is
                        specparam_winsize = 256; % scipy default
                        specparam_ovl = round(0.5*specparam_winsize);
                        specparam_nfft = 1024;
                    case 'baumgartner'
                        specparam_winfcn = @hann;
                        specparam_winsize = round(fs*(640/2048));
                        specparam_ovl = round(0.8*specparam_winsize);
                        specparam_nfft = specparam_winsize;
                    otherwise
                        error('Bad spectrum type')
                end
            else
                % read parameters from struct
                specparam_winfcn = params.wintype;
                specparam_winsize = params.winsize;
                specparam_ovl = params.ovl;
                specparam_nfft = params.nfft;
            end
            specparam_winvec = feval(specparam_winfcn, specparam_winsize);

            % compute spectrogram
            %%% using the default 'spectrumType' argument, which returns
            %%% the PSD
            [~, spec_f, spec_t, spec_psd] = spectrogram(x, specparam_winvec, specparam_ovl, specparam_nfft, fs);
            
            % determine time and frequency step
            delta_t = (specparam_winsize - specparam_ovl)/fs;
            delta_f = fs/specparam_nfft;
            
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
    
    
    %% SET/GET METHODS ====================================================
    methods
        
        % set -------------------------------------------------------------
        
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
        
        
        % get -------------------------------------------------------------
        
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
    
    
    %% METHODS - PUBLIC ===================================================
    methods (Access = public)
        
        % smooth ----------------------------------------------------------
        function obj = smooth(obj)
        % smooth spectrogram with 2D Gaussian kernel function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M = [1,2,1; 2,4,2; 1,2,1];
            psd_smooth = conv2(obj.psd, M);
            psd_smooth = psd_smooth(2:(end-1), 2:(end-1));
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
        %   The rest are identical to the "truncate" method in this class
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
            trunc_args = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
        
            % make plot
            obj.surf(ax, trunc_args{:});
            
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
            %drawnow % required, otherwise weird things happen when calling 'axis' in succesion multiple times
            %axis(ax, 'auto x');
            
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
        %   Identical to the "truncate" method in this class
        %
        % OUTPUT ARGUMENTS:
        %   ...............................................................
        %   "s" - Surface object
        %   ...............................................................
        %   
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            nargoutchk(0, 1)
            
            % check if first input is an Axes object
            n_args = numel(varargin);
            if n_args > 0 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                ax = varargin{1};
                params = varargin(2:end);
            else
                ax = gca();
                params = {};
            end
            
            % truncate spectrogram if needed
            obj_plot = obj.truncate(params{:});
            
            % prepare data for plotting
            %%% According to the MATLAB doc:
            %%%     "The time values in <t> correspond to the midpoint of
            %%%     each segment."
            %%% However, the surf() function does not plot bins - its data
            %%% points correspond to edges. Therefore, create new
            %%% variables that correspond to the start of each segment
            %%% for plotting instead. Furthermore, the last points should
            %%% be replicated to ensure that the final bins are plotted.
            if isnumeric(obj.t)
                delta_t = obj_plot.dt;
            else
                delta_t = seconds(obj_plot.dt);
            end
            
            t_plot = obj_plot.t - delta_t/2;
            t_plot = [t_plot, t_plot(end) + delta_t];
            
            f_plot = obj_plot.f - obj_plot.df/2;
            f_plot = [f_plot; f_plot(end) + obj_plot.df];
            
            psd_plot = obj_plot.psd_db;
            psd_plot = [psd_plot, psd_plot(:,end)];
            psd_plot = [psd_plot; psd_plot(end,:)];
            
            % create surface
            s = surf(ax, t_plot, f_plot, psd_plot, 'EdgeColor','none');
            
            % return surface object if needed
            if nargout > 0
                varargout = {s};
            end
        end
        
    end
end