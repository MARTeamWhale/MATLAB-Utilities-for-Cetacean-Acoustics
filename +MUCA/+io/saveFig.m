function saveFig(fig, outPath, figDims, units, varargin)
%
% Function to simplify saving figures in MATLAB.
% Successor to "saveFigInPixels".
%
% SYNTAX:
%   saveFig(fig, outPath, figDims, units)
%   saveFig(__, Name,Value)
%
% INPUT ARGUMENTS:
%   Required
%   .......................................................................
%   "fig" - handle of figure object
%   .......................................................................
%   "outPath" - path of the output file. May or may not contain the file
%       extension.
%   .......................................................................
%   "figDims" - 2-element vector containing the width and height of the
%       printed figure
%   .......................................................................
%   "units" - char string specifying the units of the figure dimensions.
%       Options include:
%           'pixels'
%           'centimeters' or 'cm'
%           'inches'
%           'points'
%           'normalized'
%   .......................................................................
%   
%   Name-Value pairs
%   .......................................................................
%   "Resolution" - integer specifying the reolution in dots per inch (DPI).
%       Default is 150. Use 300 or more for high quality.
%   .......................................................................
%   "Format" - char string specifying the output file format. If not
%       specified, will attempt to determine an appropriate format based on
%       the file extension; if a file extension was not provided, then the
%       PNG format will be used.
%   .......................................................................
%   "ForceOverwrite" - Boolean value (true/false) that determines whether
%       or not an interactive warning prompt will appear if the file being 
%       saved already exists
%   .......................................................................
%
% OUTPUT FILES:
%   image file ("outPath")
%
% DEPENDENCIES:
%   <none>
%
%
% Written by Wilfried Beslin
% Last Updated 2023-12-01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    p = inputParser;
    
    p.addRequired('fig', @(v) validateattributes(v, {'matlab.ui.Figure'}, {'scalar'}));
    p.addRequired('figPath', @(v) validateattributes(v, {'char'}, {'row'}));
    p.addRequired('figDims', @(v) validateattributes(v, {'numeric'}, {'positive','numel',2}));
    p.addRequired('units', @(v) validateattributes(v, {'char'}, {'row'}));
    p.addParameter('Resolution', 150, @(v) validateattributes(v, {'numeric'}, {'nonnegative','integer','scalar'}));
    p.addParameter('Format', '', @(v) validateattributes(v, {'char'}, {'row'}));
    p.addParameter('ForceOverwrite', false, @(v) validateattributes(v, {'logical'}, {'scalar'}));
    
    % 1) process input
    p.parse(fig, outPath, figDims, units, varargin{:});
    figDims = reshape(figDims, 1, 2);
    res = p.Results.Resolution;
    fileFormatInput = p.Results.Format;
    forceOverwrite = p.Results.ForceOverwrite;
    
    % 2) determine file format
    [outPathAdjusted, fileFormatAdjusted] = processFileFormat(outPath, fileFormatInput);
    
    % 3) prepare print arguments
    formatArg = ['-d', lower(fileFormatAdjusted)];
    resArg = ['-r', num2str(res)];
    
    % 4) check if file exists
    fileExists = isfile(outPathAdjusted);
    if fileExists && ~forceOverwrite
        questPrompt = sprintf('The file "%s" already exists.\nAre you sure you want to overwrite it?', outPathAdjusted);
        questAns = questdlg(questPrompt, 'File Exists', 'Yes', 'No', 'No');
        if ~strcmp(questAns, 'Yes')
            return
        end
    end

    % 5) get relevant figure properties
    oldProps = struct(...
        'PaperUnits', fig.PaperUnits,...
        'PaperPosition', fig.PaperPosition,...
        'PaperPositionMode', fig.PaperPositionMode);
    
    % 6) change figure properties
    figPos = [0, 0];
    if strcmpi(units, 'pixels')
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [figPos, figDims/res];
    elseif strcmpi(units, 'cm')
        fig.PaperUnits = 'centimeters';
        fig.PaperPosition = [figPos, figDims];
    else
        fig.PaperUnits = units;
        fig.PaperPosition = [figPos, figDims];
    end
    fig.PaperPositionMode = 'manual';
    
    % 7) print figure
    if fileExists
        delete(outPathAdjusted);
    end
    print(fig, outPathAdjusted, formatArg, resArg)
    
    % 8) return properties back to old
    fig.PaperUnits = oldProps.PaperUnits;
    fig.PaperPosition = oldProps.PaperPosition;
    fig.PaperPositionMode = oldProps.PaperPositionMode;
end


%% processFileFormat ------------------------------------------------------
function [figPathAdjusted, fileFormatAdjusted] = processFileFormat(figPathInput, fileFormatInput)

    % create table of supported formats
    formatCell = {...
        'jpeg',     'jpg';...
        'png',      'png';...
        'tiff',     'tif';...
        'tiffn',    'tif';...
        'meta',     'emf';...
        'pdf',      'pdf';...
        'eps',      'eps';...
        'epsc',     'eps';...
        'eps2',     'eps';...
        'epsc2',    'eps';...
        'svg',      'svg'};
    formatTable = cell2table(formatCell, 'VariableNames', {'Format','Extension'});
    extList = unique(formatTable.Extension);
    
    % process input file format
    haveFileFormat = ~isempty(fileFormatInput);
    if haveFileFormat
        % validate the input format
        fileFormatAdjusted = validatestring(fileFormatInput, formatTable.Format);
        
        % check if the file extension was added to the file path;
        % if not, then add it
        fileExt = formatTable.Extension{strcmp(fileFormatAdjusted, formatTable.Format)};
        extExpr = ['(?<=\.)(', fileExt, ')$'];
        [extMatch, iExtMatch] = regexpi(figPathInput, extExpr, 'match', 'start', 'once');
        if isempty(extMatch)
            figPathAdjusted = [figPathInput, '.', fileExt];
        else
            figPathAdjusted = [figPathInput(1:(iExtMatch-1)), lower(figPathInput(iExtMatch:end))];
        end
    else
        % look for a file extension
        extExpr = ['(?<=\.)(', strjoin(extList,'|'), ')$'];
        [extMatch, iExtMatch] = regexpi(figPathInput, extExpr, 'match', 'start', 'once');
        if isempty(extMatch)
            % if no valid extension is present, use the PNG format by default
            fileExt = 'png';
            figPathAdjusted = [figPathInput, '.', fileExt];
        else
            % use the existing extension
            fileExt = lower(extMatch);
            figPathAdjusted = [figPathInput(1:(iExtMatch-1)), lower(figPathInput(iExtMatch:end))];
        end
        
        % determine appropriate format from file extension
        % (use first instance for duplicate matches)
        iFormat = find(strcmp(fileExt, formatTable.Extension), 1, 'first');
        fileFormatAdjusted = formatTable.Format{iFormat};
    end
end