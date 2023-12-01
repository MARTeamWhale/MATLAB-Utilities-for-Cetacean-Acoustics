classdef DetectorPerformance < handle
%
% Class representing call detector performance results, given autodetection
% and ground-truth call times/frequencies, similarity index, threshold, and
% other criteria.
%
% CONSTRUCTOR SYNTAX:
%   obj = DetectorPerformance(det_times, call_times, f_max, sim_idx, sim_th)
%   obj = DetectorPerformance(det_times, det_freqs, call_times, call_freqs, sim_idx, sim_th)
%
% CONSTRUCTOR INPUT ARGUMENTS:
%   .......................................................................
%   "det_times" - nDet-by-2 matrix representing the start and end times of
%       each detection
%   .......................................................................
%   "det_freqs" - nDet-by-2 matrix representing the minimum and maximum
%       frequency of each detection. If not specified, frequencies will
%       range from 0 to f_max for all detections.
%   .......................................................................
%   "call_times" - nCalls-by-2 matrix representing the start and end times
%       of each true call in the dataset
%   .......................................................................
%   " call_freqs" - nCalls-by-2 matrix representing the minimum and maximum
%       frequency of each true call. If not specified, frequencies will
%       range from 0 to f_max for all calls.
%   .......................................................................
%   "f_max" - maximum frequency. Specify this only if frequency ranges for
%       individual detections and calls are not available or not important.
%   .......................................................................
%   "sim_idx" - char string representing the similarity index to use when
%       comparing detections against true calls. Options are:
%           'ovl', 'jaccard', or 'dice'.
%   .......................................................................
%   "sim_th" - similarity threshold above which a call-detection pair must 
%       score to be considered a match. Must be a number between 0 and 1.
%   .......................................................................
%
% DEPENDENCIES:
%   MUCA.dcs_analysis.ConfusionMatrix ["makeConfusionMatrix" method]
%
%
% Written by Wilfried Beslin
% Last updated 2023-11-30
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% - The similarity matrices use a rudimentary sparse matrix implementation
% for better memory efficiency, but it could likely be improved further

    properties
        sim_idx         % similarity index; may be either 'ovl', 'jaccard', or 'dice'
        sim_th          % similarity threshold; value between 0 and 1
    end
    properties (SetAccess = private)
        detections      % table of detections with fields 'StartTime', 'EndTime', 'MinFreq', 'MaxFreq'
        calls           % table of ground-truth calls with fields 'StartTime', 'EndTime', 'MinFreq', 'MaxFreq'
        excluded_calls  % vector of logicals (one for each call) indicating if a call is to be excluded from the final score or not (i.e., matches will not be counted as either positive or negative)
    end
    properties (SetAccess = private, Dependent)
        num_detections
        num_calls
        similarity              % N-by-M matrix of similarity scores, where N is number of detections and M is number of true calls. Depends on sim_idx.
        call_match              % vector indicating which call each detection is matched to, if any. Depends on similarity and sim_th.
        call_match_adjusted     % like call_match, but only one detection per call. In the event there are multiple dets per call, the det with highest similarity wins, and the others are marked as NaN.
        scores                  % 0 if no match, 1 if match, NaN if redundant or matching an excluded call. Depends on call_match_adjusted and call_match.
        missed_calls            % logical indices of calls for which there were no matches
        broken_calls            % logical indices of calls for which there are multiple matches
    end
    properties (Constant, Hidden)
        sim_fcn_list = struct(...
            'ovl', @WSL.detection.DetectorPerformance.calculate_ovl,...
            'jaccard', @WSL.detection.DetectorPerformance.calculate_jaccard,...
            'dice', @WSL.detection.DetectorPerformance.calculate_dice);
        sim_idx_list = fieldnames(WSL.detection.DetectorPerformance.sim_fcn_list);
    end
    properties (SetAccess = private, Hidden)
        similarity_all  % struct of sparse similarity matrices for each similarity index
    end
    
    % Constructor ---------------------------------------------------------
    methods
        function self = DetectorPerformance(varargin)
        
            if nargin == 5
                [det_times, call_times, f_max, sim_idx, sim_th] = varargin{:};
                det_freqs = [];
                call_freqs = [];
            elseif nargin == 6
                [det_times, det_freqs, call_times, call_freqs, sim_idx, sim_th] = varargin{:};
            else
                error('Incorrect number of input arguments')
            end
            
            % validation
            assert(size(det_times,2) == 2 && size(call_times,2) == 2, 'Time intervals must be N-by-2 matrices')
            assert(all(det_times(:,1) <= det_times(:,2)) && all(call_times(:,1) <= call_times(:,2)), 'Start times cannot be greater than end times')
            if isempty(det_freqs)
                det_freqs = repmat([0,f_max], size(det_times,1), 1);
                call_freqs = repmat([0,f_max], size(call_times,1), 1);
            else
                assert(all(size(det_freqs) == size(det_times)) && all(size(call_freqs) == size(call_times)), 'Time and frequency matrices must have identical size')
                assert(all(det_freqs(:,1) <= det_freqs(:,2)) && all(call_freqs(:,1) <= call_freqs(:,2)), 'Min frequencies cannot be greater than max frequencies')
            end
            
            % create properties
            self.detections = table(det_times(:,1), det_times(:,2), det_freqs(:,1), det_freqs(:,2), 'VariableNames',{'StartTime','EndTime','MinFreq','MaxFreq'});
            self.calls = table(call_times(:,1), call_times(:,2), call_freqs(:,1), call_freqs(:,2), 'VariableNames',{'StartTime','EndTime','MinFreq','MaxFreq'}); 
            self.sim_idx = sim_idx;
            self.sim_th = sim_th;
            self.excluded_calls = false(size(self.calls));
            
            self.computeSimilarity();
        end
    end
    
    %% SET METHODS ========================================================
    methods
        
        % set.sim_idx -----------------------------------------------------
        function set.sim_idx(self, v)
            self.sim_idx = validatestring(v, self.sim_idx_list);
        end
        
        % set.sim_th ------------------------------------------------------
        function set.sim_th(self, v)
            validateattributes(v,{'numeric'},{'scalar','>=',0,'<=',1})
            self.sim_th = v;
        end
        
    end
    
    %% GET METHODS ========================================================
    methods
        
        % get.num_detections ----------------------------------------------
        function v = get.num_detections(self)
            v = height(self.detections);
        end
        
        % get.num_calls ---------------------------------------------------
        function v = get.num_calls(self)
            v = height(self.calls);
        end
        
        % get.similarity --------------------------------------------------
        function v = get.similarity(self)
            v = full(self.similarity_all.(self.sim_idx));
        end
        
        % get.call_match --------------------------------------------------
        function v = get.call_match(self)
            [sim_max, i_max] = max(self.similarity, [], 2, 'omitnan');
            i_max(sim_max < self.sim_th | isnan(sim_max)) = NaN;
            v = i_max;
        end
        
        % get.call_match_adjusted -----------------------------------------
        function v = get.call_match_adjusted(self)
            match_raw = self.call_match;
            sim_max = max(self.similarity, [], 2, 'omitnan');
            matched_calls = unique(match_raw(~isnan(match_raw)))'; % list of calls with a match
            v = NaN(size(match_raw));
            for ii = matched_calls
                match_ii = match_raw == ii;  % logical index of detections matching current call
                sim_max_ii = sim_max.*match_ii;  % vector of max similarities, with all non-target detections set to zero
                [abs_sim_max_ii, i_abs_sim_max_ii] = max(sim_max_ii, [], 'omitnan');
                v(i_abs_sim_max_ii) = abs_sim_max_ii;
            end
        end
        
        % get.scores ------------------------------------------------------
        function v = get.scores(self)
            raw_scores = ~isnan(self.call_match);
            v = double(raw_scores);
            
            % turn redundant matches into NaN
            is_redundant = raw_scores & isnan(self.call_match_adjusted);
            v(is_redundant) = NaN;
            
            % turn matches with excluded calls into NaNs
            matching_excluded = ismember(self.call_match, find(self.excluded_calls));
            v(matching_excluded) = NaN;
        end
        
        % get.missed_calls ------------------------------------------------
        function v = get.missed_calls(self)
            i_calls = (1:self.num_calls)';
            v = ~ismember(i_calls, self.call_match);
        end
        
        % get.broken_calls ------------------------------------------------
        function v = get.broken_calls(self)
            i_calls = (1:self.num_calls)';
            
            match_raw = self.call_match;
            match_adjusted = self.call_match_adjusted;
            is_broken_part = ~isnan(match_raw) & isnan(match_adjusted);
            i_broken = unique(match_raw(is_broken_part));
            
            v = ismember(i_calls, i_broken);
        end
        
    end
    
    %% PUBLIC METHODS =====================================================
    methods
        
        % assignExcludedCalls ---------------------------------------------
        function assignExcludedCalls(self, exclude)
        % Set which calls are to be excluded, if any
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            validateattributes(exclude,{'logical'},{'size',[self.num_calls,1]})
            self.excluded_calls = exclude;
        end
        
        % makeConfusionMatrix ---------------------------------------------
        function cm = makeConfusionMatrix(self)
        % Produces a ConfusionMatrix object from current properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            import MUCA.dcs_analysis.ConfusionMatrix
            
            % actual
            a = self.scores;  
            a(isnan(a)) = [];  % remember, NaNs are redundant detections and detections with excluded calls
            a = logical(a);
            
            % predicted
            p = true(size(a));
            
            % add missed calls at end
            num_missed = sum(self.missed_calls(~self.excluded_calls));
            a = [a; true(num_missed,1)];
            p = [p; false(num_missed,1)];
            
            % generate matrix object
            cm = ConfusionMatrix(a,p);
        end
        
        % plot ------------------------------------------------------------
        function varargout = plot(self, varargin)
        % Plot the autodetection and ground-truth polygons
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            nargoutchk(0,1)
            narginchk(1,2)
            
            if isempty(varargin)
                ax = gca;
                cla(ax)
            else
                ax = varargin{1};
            end
            ax.NextPlot = 'add';
            
            x_det_plot = [self.detections.StartTime, self.detections.StartTime, self.detections.EndTime, self.detections.EndTime]';
            y_det_plot = [self.detections.MinFreq, self.detections.MaxFreq, self.detections.MaxFreq, self.detections.MinFreq]';

            x_calls_plot = [self.calls.StartTime, self.calls.StartTime, self.calls.EndTime, self.calls.EndTime]';
            y_calls_plot = [self.calls.MinFreq, self.calls.MaxFreq, self.calls.MaxFreq, self.calls.MinFreq]';

            col = lines(2);

            hpa = fill(ax, x_calls_plot, y_calls_plot, col(1,:), 'DisplayName','True Calls');
            hpb = fill(ax, x_det_plot, y_det_plot, col(2,:), 'FaceAlpha',0.5, 'DisplayName','Detections');
            legend(ax, [hpa(1),hpb(1)])
            xlabel(ax, 'Time')
            ylabel(ax, 'Frequency')
            grid(ax, 'on');
            box(ax, 'on');
            
            % return axes if specified
            if nargout == 1
                varargout{1} = ax;
            end
        end
    end
    
    
    %% STATIC METHODS =====================================================
    methods (Static)
        
        % calculate_ovl ---------------------------------------------------
        function ovl = calculate_ovl(xa, ya, xb, yb)
        % Calculate overlap coefficient between pairs of 2D square volumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            vnab = WSL.detection.DetectorPerformance.intersection2(xa, ya, xb, yb);
            va = WSL.detection.DetectorPerformance.area2(xa, ya);
            vb = WSL.detection.DetectorPerformance.area2(xb, yb);
            
            ovl = vnab./min([va,vb],[],2);
        end
        
        % calculate_jaccard -----------------------------------------------
        function jac = calculate_jaccard(xa, ya, xb, yb)
        % Calculate Jaccard similarity indices between pairs of 2D square 
        % volumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            vnab = WSL.detection.DetectorPerformance.intersection2(xa, ya, xb, yb);
            vuab = WSL.detection.DetectorPerformance.union2(xa, ya, xb, yb);
            
            jac = vnab./vuab;
        end
        
        % calculate_dice --------------------------------------------------
        function dsc = calculate_dice(xa, ya, xb, yb)
        % Calculate Sorensen-Dice similarity coefficients between pairs of
        % 2D square volumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            vnab = WSL.detection.DetectorPerformance.intersection2(xa, ya, xb, yb);
            va = WSL.detection.DetectorPerformance.area2(xa, ya);
            vb = WSL.detection.DetectorPerformance.area2(xb, yb);

            dsc = (2*vnab)./(va + vb);
        end
        
        % area2 -----------------------------------------------------------
        function a = area2(x,y)
        % Calculate areas of 2D square volumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            a = diff(x,[],2).*diff(y,[],2);
        end
        
        % intersection2 ---------------------------------------------------
        function n = intersection2(xa, ya, xb, yb)
        % calculate intersections between pairs of 2D square volumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            nx = min([xa(:,2),xb(:,2)],[],2) - max([xa(:,1),xb(:,1)],[],2);
            nx(nx < 0) = 0;
            ny = min([ya(:,2),yb(:,2)],[],2) - max([ya(:,1),yb(:,1)],[],2);
            ny(ny < 0) = 0;
            n = nx.*ny;
        end
        
        % union2 ----------------------------------------------------------
        function u = union2(xa, ya, xb, yb)
        % Calculate union between pairs of 2D square volumes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            u = WSL.detection.DetectorPerformance.area2(xa,ya) + WSL.detection.DetectorPerformance.area2(xb,yb) - WSL.detection.DetectorPerformance.intersection2(xa,ya,xb,yb);
        end
        
    end
    
    
    %% PRIVATE METHODS ====================================================
    methods (Access = private)
        
        % computeSimilarity -----------------------------------------------
        function computeSimilarity(self)
        % Calculate the similarity matrix between autodetection and true
        % call time-frequency boxes, for all similarity indices.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            disp('DetectorPerformance: calculating similarity matrix')
        
            n = self.num_detections;
            m = self.num_calls;
            k = numel(self.sim_idx_list);
            det_times = [self.detections.StartTime, self.detections.EndTime];
            det_freqs = [self.detections.MinFreq, self.detections.MaxFreq];
            call_times = [self.calls.StartTime, self.calls.EndTime];
            call_freqs = [self.calls.MinFreq, self.calls.MaxFreq];
            %sim_fcn = self.sim_fcn_list.(self.sim_idx);
            
            % if times are datetimes, convert to doubles - more efficient
            if isdatetime(det_times)
                dt_ref = min([det_times(1),call_times(1)]);
                det_times = seconds(det_times - dt_ref);
                call_times = seconds(call_times - dt_ref);
            end
            
            % initialize container
            s = struct();
            
            % loop through each index
            for ii = 1:k
                sim_idx_ii = self.sim_idx_list{ii};
                sim_fcn_ii = self.sim_fcn_list.(sim_idx_ii);
                
                s_ii = zeros(n,m, 'single');
                
                % loop through each detection and compare against calls
                for jj = 1:n
                    fprintf('    index %d (%s): detection %d/%d\n', ii, sim_idx_ii, jj, n)
                    det_time_array_jj = repelem(det_times(jj,:), m, 1);
                    det_freq_array_jj = repelem(det_freqs(jj,:), m, 1);
                    s_ii(jj,:) = feval(sim_fcn_ii, det_time_array_jj, det_freq_array_jj, call_times, call_freqs)';
                end
                
                % convert to sparse matrix
                s.(sim_idx_ii) = sparse(double(s_ii));
            end
            
            % assign results to object property
            self.similarity_all = s;
        end
   
    end
end