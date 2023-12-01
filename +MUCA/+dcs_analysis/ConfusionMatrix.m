classdef ConfusionMatrix
%
% Value class representing a binary confusion matrix. Input consists of
% the outcome of several yes/no trials, for both the actual case and the
% predicted case.
%
% CONSTRUCTOR SYNTAX:
%   obj = ConfusionMatrix(actualClass, predictedClass)
%   obj = ConfusionMatrix(__, classNames)
%
% CONSTRUCTOR INPUT ARGUMENTS:
%   .......................................................................
%   "actualClass" - vector of true/false values indicating the actual class
%       of each instance
%   .......................................................................
%   "predictedClass" - vector of true/false values indicating the predicted
%       class of each instance
%   .......................................................................
%   "classNames" - 2-element cell array of strings corresponding the names
%       of each class: The first element is class 1, and the second is
%       class 0. Default is {'positive', 'negative'}.
%   .......................................................................
%
% DEPENDENCIES:
%   Statistics and Machine Learning Toolbox ["confusionchart" method]
%   Neural Network Toolbox ["plotconfusion" method]
%
%
% Written by Wilfried Beslin Jun. 2018
% Last updated 2023-11-30
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    properties (SetAccess = private)
        caseNames   % 2-element cellstr; corresponds to the names of positives and negatives
    end
    properties (SetAccess = private, Dependent)
        labels % table of actual and predicted class labels (based on caseNames)
        nTrials % total number of data points
        P   % number of condition positives
        N   % number of condition negatives
        TP  % number of true positives: those predicted positives that are actually positive
        TN  % number of true negatives: those predicted negatives that are actually negative
        FP  % number of false positives: those predicted positives that are actually negative
        FN  % number of false negatives: those predicted negatives that are actually positive
        TPR % true positive rate: proportion of condition positives that were predicted positive
        TNR % true negative rate: proportion of condition negatives that were predicted negative
        FPR % false positive rate: proportion of condition negatives that were predicted positive
        FNR % false negative rate: proportion of condition positives that were predicted negative
        PPV % positive predictive value: proportion of predicted positives that are actually positive
        NPV % negative predictive value: proportion of predicted negatives that are actually negative
        FDR % false discovery rate: proportion of predicted positives that are actually negative
        FOR % false ommision rate: proportion of predicted negatives that are actually positive
        accuracy % proportion of correct predictions
        waccuracy % weighted accuracy; average of sensitivity and specificity
        sensitivity % alias for TPR
        specificity % alias for TNR
        precision % alias for PPV
        recall % alias for TPR
        F1Score % harmonic mean of precision and recall
    end
    properties (SetAccess = private, Hidden)
        data    % table of yes/no trial results for actual and predicted case
    end
    
    % Constructor ---------------------------------------------------------
    methods
        function obj = ConfusionMatrix(varargin)
            if nargin > 0
                % parse input
                p = inputParser;
                validTrial = @(arg) validateattributes(arg,{'logical'},{'vector'});
                %%% actual
                p.addRequired('actual',validTrial)
                %%% predicted
                p.addRequired('predicted',validTrial)
                %%% names
                validNames = @(arg) iscellstr(arg) && numel(arg) == 2;
                defaultNames = {'positive','negative'};
                p.addOptional('names',defaultNames,validNames)
                %%% parsing
                p.parse(varargin{:})
                actual = p.Results.actual;
                predicted = p.Results.predicted;
                names = p.Results.names;
                % end input parsing

                % reshape vectors
                actual = reshape(actual,numel(actual),1);
                predicted = reshape(predicted,numel(predicted),1);
                names = reshape(names,1,2);

                % create table
                dataTable = table(actual,predicted,'VariableNames',{'Actual','Predicted'});

                % assign properties
                obj.data = dataTable;
                obj.caseNames = names;
            end
        end
    end
    
    %% GET METHODS ========================================================
    methods
        % get.labels ------------------------------------------------------
        function c = get.labels(obj)
            cCell = repmat({''}, obj.nTrials, 2);
            dMat = obj.data{:,:};
            
            cCell(dMat) = obj.caseNames(1);
            cCell(~dMat) = obj.caseNames(2);
            
            c = cell2table(cCell, 'VariableNames',obj.data.Properties.VariableNames);
        end
        
        % get.nTrials -----------------------------------------------------
        function n = get.nTrials(obj)
            n = height(obj.data);
        end
        
        % get.P -----------------------------------------------------------
        function c = get.P(obj)
            c = sum(obj.data.Actual);
        end
        
        % get.N -----------------------------------------------------------
        function c = get.N(obj)
            c = sum(~obj.data.Actual);
        end
        
        % get.TP ----------------------------------------------------------
        function c = get.TP(obj)
            c = sum(obj.data.Predicted & obj.data.Actual);
        end
        
        % get.TN ----------------------------------------------------------
        function c = get.TN(obj)
            c = sum(~obj.data.Predicted & ~obj.data.Actual);
        end
        
        % get.FP ----------------------------------------------------------
        function c = get.FP(obj)
            c = sum(obj.data.Predicted & ~obj.data.Actual);
        end
        
        % get.FN ----------------------------------------------------------
        function c = get.FN(obj)
            c = sum(~obj.data.Predicted & obj.data.Actual);
        end
        
        % get.TPR ---------------------------------------------------------
        function r = get.TPR(obj)
            r = obj.TP/obj.P;
        end
        
        % get.TNR ---------------------------------------------------------
        function r = get.TNR(obj)
            r = obj.TN/obj.N;
        end
        
        % get.FPR ---------------------------------------------------------
        function r = get.FPR(obj)
            r = obj.FP/obj.N;
        end
        
        % get.FNR ---------------------------------------------------------
        function r = get.FNR(obj)
            r = obj.FN/obj.P;
        end
        
        % get.PPV ---------------------------------------------------------
        function r = get.PPV(obj)
            r = obj.TP/(obj.TP + obj.FP);
        end
        
        % get.NPV ---------------------------------------------------------
        function r = get.NPV(obj)
            r = obj.TN/(obj.TN + obj.FN);
        end
        
        % get.FDR ---------------------------------------------------------
        function r = get.FDR(obj)
            r = obj.FP/(obj.FP + obj.TP);
        end
        
        % get.FOR ---------------------------------------------------------
        function r = get.FOR(obj)
            r = obj.FN/(obj.FN + obj.TN);
        end
        
        % get.accuracy ----------------------------------------------------
        function a = get.accuracy(obj)
            a = (obj.TP + obj.TN)/obj.nTrials;
        end
        
        % get.waccuracy ---------------------------------------------------
        function wa = get.waccuracy(obj)
            wa = mean([obj.sensitivity,obj.specificity]);
        end
        
        % get.sensitivity -------------------------------------------------
        function r = get.sensitivity(obj)
            r = obj.TPR;
        end
        
        % get.specificity -------------------------------------------------
        function r = get.specificity(obj)
            r = obj.TNR;
        end
        
        % get.precision ---------------------------------------------------
        function r = get.precision(obj)
            r = obj.PPV;
        end
        
        % get.recall ------------------------------------------------------
        function r = get.recall(obj)
            r = obj.TPR;
        end
        
        % get.F1Score -----------------------------------------------------
        function s = get.F1Score(obj)
            s = (2.*obj.precision.*obj.recall)./(obj.precision + obj.recall);
        end
    end
    
    
    %% PUBLIC METHODS =====================================================
    methods
        % confusionchart --------------------------------------------------
        function cm = confusionchart(obj, varargin)
        % Overload of Statistics and Machine Learning Toolbox function.
        % Works just like the original, except data does not need to be
        % specified (use the object instead).
        % 
        % Recommend using the following parameters:
        % - 'RowSummary','row-normalized'
        % - 'ColumnSummary','column-normalized'
        %
        % Using these settings, you end up with the following:
        % - Right summary matrix:
        %   [ TNR,              FPR;
        %     TPR (RECALL),     FNR (MISS RATE) ]
        % Bottom summary matrix:
        %   [ Negative Predictive Value,    Positive Predictive Value (PRECISION);
        %     False Omission Rate,          False Discovery Rate ]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            lb = obj.labels;
            if ~isempty(varargin) && isgraphics(all(varargin{1}))
                cm = confusionchart(varargin{1}, lb.Actual, lb.Predicted, varargin{2:end});
            else
                cm = confusionchart(lb.Actual, lb.Predicted, varargin{:});
            end
        end
        
        % plotconfusion ---------------------------------------------------
        function varargout = plotconfusion(obj,varargin)
        % Plots the confusion matrix.
        % Overload of Neural Network Toolbox function.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % validate output
            nargoutchk(0,1)
            
            % parse input
            if nargin > 0 && all(isgraphics(varargin{1}))
                hAx = varargin{1};
                iMainArgs = true(1,numel(varargin));
                iMainArgs(1) = false;
                mainArgs = varargin(iMainArgs);
            else
                hAx = [];
                mainArgs = varargin;
            end
            p = inputParser;
            %%% name
            validName = @(arg) ischar(arg) || iscellstr(arg);
            defaultName = '';
            p.addOptional('name',defaultName,validName)
            %%% parse
            p.parse(mainArgs{:})
            name = p.Results.name;
            % end input parsing
            
            % "plotconfusion" takes 3 input arguments: "target", "output",
            % and "name". "target" and "output" are binary matrices where
            % rows = classes, columns = observations. "name" is an optional 
            % string (or cellstr) add to the title of the plot.
            %%% target
            targetPresent = obj.data.Actual';
            targetAbsent = ~targetPresent;
            target = [targetPresent;targetAbsent];
            %%% output
            outputPresent = obj.data.Predicted';
            outputAbsent = ~outputPresent;
            output = [outputPresent;outputAbsent];
            %%% classes
            classLabels = [obj.caseNames,{''}];
            
            % create plot
            if isempty(hAx)
                % create new figure
                figure();
            else
                % set current axes and figure
                hFig = hAx.Parent;
                set(groot,'CurrentFigure',hFig);
                set(hFig,'CurrentAxes',hAx);
            end
            plotconfusion(target,output,name)
            hAx = gca;
            hAx.XTickLabel = classLabels;
            hAx.YTickLabel = classLabels;
            
            % set output if requested
            if nargout > 0
                varargout{1} = hAx;
            end
        end
    end
end