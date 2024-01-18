function dt = readDateTime(fileNames, varargin)
%
% Returns a datetime object corresponding to the timestamp in a file name
% (or list of file names). Full paths are also accepted. Note that in case
% multiple timestamps are found in the path, only the last match is
% returned.
%
% IMPORTANT: this only works if the date/time is encoded somewhere in the
% filename as "yyyyMMdd.HHmmss" or "yyyyMMdd.HHmmss.SSS", where "." is any
% character that is not a digit or a slash, and SSS represents
% milliseconds. The presence or absence of milliseconds is determined
% automatically.
%
% SYNTAX:
%   dt = readDateTime(fileNames)
%   dt = readDateTime(fileNames, dtFormat)
%   dt = readDateTime(__, 'SuppressWarnings',value)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "fileNames" - char string or cell array of char strings containing the
%       names of files with timestamps
%   .......................................................................
%   "dtFormat" [OPTIONAL] - char string representing the date-time display 
%       format to use (see "https://www.mathworks.com/help/matlab/ref/datetime.html#buhzxmk-1-Format").
%       If not specified, the format will be determined automatically using
%       a simplified ISO 8601 format that may or may not display
%       milliseconds, depending on whether this information is included in
%       the timestamps or not.
%   .......................................................................
%   "SuppressWarnings" [OPTIONAL, NAME-VALUE PAIR] - logical value that
%       determines if warnings should be hidden or not. Default is false. 
%       Warnings appear if timestamps could not be read from filenames.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "dt" - array of datetime objects corresponding to the file timestamps
%   .......................................................................
%
% OUTPUT FILES:
%   <none>
%
% DEPENDENCIES:
%   MUCA.time.isoFormat
%
%
% Written by Wilfried Beslin
% Last Updated 2024-01-18 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% 2022-05-03
% - there could be an argument to decide which date-times are
% returned, if multiple matches are found
% -- However, this proved more work than expected to implement
    
    import MUCA.time.isoFormat
    
    % parse input
    p = inputParser();
    p.addRequired('fileNames', @(val)ischar(val)||iscellstr(val))
    if nargin == 2 || nargin == 4
        p.addRequired('dtFormat', @ischar)
    end
    p.addParameter('SuppressWarnings', false, @(val)validateattributes(val,{'logical'},{'scalar'}))
    
    p.parse(fileNames,varargin{:});
    if ismember('dtFormat',fieldnames(p.Results))
        dtFormat = p.Results.dtFormat;
    else
        dtFormat = '';
    end
    showWarnings = ~p.Results.SuppressWarnings;
    
    % convert fileNames to cellstr if it isn't already
    fileNames = cellstr(fileNames);
    numFiles = numel(fileNames);
    
    % define regex for reading datetime
    %%% This expression finds strings of the form "DDDDDDDD.DDDDDD(.DDD)", 
    %%% where the stuff in brackets is optional and "." is any character 
    %%% that is not a digit or a slash
    expr = '(?<!\d)\d{8}[^\d\\/]\d{6}([^\d]\d{3})?(?!\d)';
    
    % match datetime strings in filenames
    dtStrs = regexp(fileNames, expr, 'match');
    
    % process each file
    dt = NaT(size(fileNames)); % initialize datetime array
    has_ms = false(numFiles, 1);
    for ii = 1:numFiles
        % get all matches for current file
        dtStrs_ii = [dtStrs{ii}];  
        
        try
            % isolate only the last match (to ensure that only the match
            % for filename is used, in case the input is a full file path
            % that may contain other dates-times
            dtStr_ii = dtStrs_ii{end};
        
            % extract datetime from filename
            [dt_ii, has_ms_ii] = readDTStr(dtStr_ii);
            
            dt(ii) = dt_ii;
            has_ms(ii) = has_ms_ii;
        catch
            % unable to read datetime
            if showWarnings
                warning('Could not read timestamp for file "%s"', fileNames{ii})
            end
        end
    end
    
    % set format
    if isempty(dtFormat)
        isoArgs = {'simplified'};
        if any(has_ms)
            isoArgs = [isoArgs,{'milliseconds'}];
        end
        dtFormat = isoFormat(isoArgs{:});
    end
    try
        dt.Format = dtFormat;
    catch
        warning('The specified date-time display format is not supported')
    end
end


%% readDTStr --------------------------------------------------------------
function [dt, has_ms] = readDTStr(dtStr)
% Extract datetime from char string and return as datetime object.
% Assumes the strings are of the form "DDDDDDDD.DDDDDD(.DDD)"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define indices of separating characters
    iSep1 = 9;
    iSep2 = 16;
    
    % extract datetime info
    if length(dtStr) < iSep2
        has_ms = false;
        dtStrSub = dtStr([1:(iSep1-1),(iSep1+1):end]);
        dtInputFormat = 'yyyyMMddHHmmss';
    else
        has_ms = true;
        dtStrSub = dtStr([1:(iSep1-1), (iSep1+1):(iSep2-1), (iSep2+1):end]);
        dtInputFormat = 'yyyyMMddHHmmssSSS';
    end
    dt = datetime(dtStrSub, 'InputFormat',dtInputFormat);
end