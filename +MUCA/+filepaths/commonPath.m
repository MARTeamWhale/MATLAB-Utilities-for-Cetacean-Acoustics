function rootPath = commonPath(filePaths)
%
% Look for a common path among a list of file paths.
%
% SYNTAX:
%   rootPath = commonPath(filePaths)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "filePaths" - cell array of file path strings
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "rootPath" - common path among all items in list
%   .......................................................................
%
% OUTPUT FILES:
%   <none>
%
% DEPENDENCIES:
%   <none>
%
%
% Written by Wilfried Beslin
% Last Updated 2022-11-23 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    filePaths_char = char(filePaths);
    numCols = size(filePaths_char, 2);
    
    % loop through each character column and check if they are the same or
    % not
    idx = 1;
    idxLastFileSep = NaN;
    allCharsEqual = true;
    while allCharsEqual && idx <= numCols
        uniqueChars = unique(filePaths_char(:,idx));
        allCharsEqual = numel(uniqueChars) == 1;
        
        % make note if the current character is a file seperator
        if allCharsEqual && uniqueChars == filesep
            idxLastFileSep = idx;
        end
        
        % increment index
        idx = idx + 1;
    end
    
    % isolate the part that is common to all paths before the last file
    % seperator, if any
    if ~isnan(idxLastFileSep)
        rootPath = filePaths_char(1, (1:(idxLastFileSep-1)));
    else
        % no common paths
        rootPath = '';
    end
end