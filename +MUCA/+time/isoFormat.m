function [ft, ftDate, ftTime] = isoFormat(varargin)
%
% Returns a char string representing a date-time display format based on 
% the ISO 8601 standard, for use with MATLAB 'datetime' objects.
%
% SYNTAX:
%   isoFormat()
%   isoFormat(option1, option2, ... optionN)
%
% OPTIONAL INPUT ARGUMENTS (char strings):
%   .......................................................................
%   'long' - use the ISO 8601 extended format (i.e, 
%       yyyy-MM-ddTHH:mm:ss). This is the default. Cannot be specified 
%       at the same time as 'compact'.
%   .......................................................................
%   'compact' - use the ISO 8601 basic format (i.e., yyyyMMddTHHmmss).
%       Cannot be specified at the same time as 'long'.
%   .......................................................................
%   'simplified' - replace the T seperator charactor with a space (for long
%       format) or an underscore (for compact format). Note that this
%       deviates a bit from the regular ISO format, but it often makes 
%       things easier to read.
%   .......................................................................
%   'milliseconds' - add fractional seconds, to three decimal places (e.g.,
%       yyyy-MM-ddTHH:mm:ss.SSS).
%   .......................................................................
%   'utc' - emphasize that a time is in UTC by adding 'Z' at the end.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   ft - char string representing the full datetime format
%   .......................................................................
%   ftDate - char string representing only the date component, e.g., 
%       yyyy-MM-dd
%   .......................................................................
%   ftTime - char string representing only the time component, e.g. 
%       HH:mm:ss. NOTE: to use this with a duration object, convert it to 
%       lowercase first.
%   .......................................................................
%
% OUTPUT FILES:
%   <none>
%
% DEPENDENCIES:
%   <none>
%
% NOTES:
% - information on the ISO 8601 standard can be found here: 
% https://en.wikipedia.org/wiki/ISO_8601 
%
%
% Written by Wilfried Beslin
% Last updated 2023-12-01, using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    assert(iscellstr(varargin), 'Input must only consist of char strings!')

    validStyles = {'long', 'compact'};
    simplifyFlag = 'simplified';
    millisecFlag = 'milliseconds';
    utcFlag = 'utc';
    
    possibleOptions = [validStyles, {simplifyFlag, millisecFlag, utcFlag}];
    
    userOptions = unique(lower(varargin));
    if numel(userOptions) ~= numel(varargin)
        warning('Ignoring duplicate entries')
    end
    
    assert(all(ismember(userOptions, possibleOptions)), 'One or more input options is not recognized.\nCheck the function''s ''help'' for a list of valid options.')
    
    % determine format style (long or compact)
    userStyle = intersect(validStyles, userOptions);
    switch numel(userStyle)
        case 0
            formatStyle = 'long';
        case 1
            formatStyle = userStyle{:};
        otherwise
            error('Only one format style can be specified!')
    end
    
    % determine if simplify flag is active
    if ismember(simplifyFlag, userOptions)
        simplify = true;
    else
        simplify = false;
    end
    
    % determine if millisecond flag is active
    if ismember(millisecFlag, userOptions)
        addMilliseconds = true;
    else
        addMilliseconds = false;
    end
    
    % determine if UTC flag is active
    if ismember(utcFlag, userOptions)
        addUTC = true;
    else
        addUTC = false;
    end
    
    % END INPUT PARSING
    
    % begin with long datetime format, and make adjustments depending on
    % options
    ft = 'yyyy-MM-dd''T''HH:mm:ss';
    
    if strcmp(formatStyle, 'compact')
        ft = erase(ft, {'-',':'});
    end
    
    if simplify
        simpleSep = struct('long',' ', 'compact','_');
        ft = strrep(ft,'''T''', simpleSep.(formatStyle));
    end
    
    if addMilliseconds
        ft = [ft,'.SSS'];
    end
    
    if addUTC
        ft = [ft,'''Z'''];
    end
    
    % get individual date and time formats
    splitExpr = '(y[yMd-]*d)|(H[HmsS\.:]*(s|S))';
    splitMatch = regexp(ft, splitExpr, 'match');
    ftDate = splitMatch{1};
    ftTime = splitMatch{2};
end