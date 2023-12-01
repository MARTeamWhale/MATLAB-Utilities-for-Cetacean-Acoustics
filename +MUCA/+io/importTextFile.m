function fileStr = importTextFile(filePath, varargin)
%
% Read a text file as a cell array of strings, where each cell corresponds
% to one row in the file.
%
% SYNTAX:
%   fileStr = importTextFile(filePath)
%   fileStr = importTextFile(filePath, trimLeadingWhitespace)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "filePath" - char string specifying path to text file
%   .......................................................................
%   "trimLeadingWhitespace" - true/false value specifying whether
%       whitespace at the beginning of each line should be removed or not.
%       Default is true.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "fileStr" - cell array of strings containing the information within the
%       input text file
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
% Last Updated 2023-12-01 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    narginchk(1, 2)
    if isempty(varargin)
        trimLeadingWhitespace = true;
    else
        trimLeadingWhitespace = varargin{1};
    end
    
    if trimLeadingWhitespace
        whitespaceArgs = {};
    else
        whitespaceArgs = {'Whitespace', ''};
    end

    fileStr = {};
    errMsg = '';

    % open file (read-only mode)
    [fileID, readErr] = fopen(filePath, 'r');
    
    if fileID ~= -1
        try
            % read file contents as text
            fileStr = textscan(fileID, '%s', 'Delimiter',{'\n','\r'}, whitespaceArgs{:});
            fileStr = fileStr{:};
        catch ME
            errMsg = sprintf('Failed to read contents of %s\n%s', filePath, ME.message);
        end
        
        % close file
        fclose(fileID);
    else
        errMsg = sprintf('Unable to open %s\n%s', filePath, readErr);
    end
    
    % issue warning if there were problems
    if isempty(fileStr)
        warning(errMsg)
    end
end