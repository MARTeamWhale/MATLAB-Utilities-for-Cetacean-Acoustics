function errMsg = writeTextFile(txt, filePath, varargin)
%
% Write a char string or cell array of strings to a text file.
%
% SYNTAX:
%   errMsg = writeTextFile(txt, filePath)
%   errMsg = writeTextFile(txt, filePath, existingFileOpt)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "txt" - char string (or cell array of char strings) containing the text
%       to write to file
%   .......................................................................
%   "filePath" - char string specifying the path of the output file
%   .......................................................................
%   "existingFileOpt" - char string specifying the action to take in case
%       the output file already exists. The options are:
%           'prompt' [DEFAULT] - prompt user for an action
%           'overwrite' - overwrite the file
%           'append' - append text to the existing file
%           'cancel' - do not write anything
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "errMsg" - char string describing any errors that were encountered when
%       trying to write to the file
%   .......................................................................
%
% OUTPUT FILES:
%   text file ("filePath")
%
% DEPENDENCIES:
%   <none>
%
% NOTES:
% - When appending to files, the new text will be added immediately after
% the last character in the end of the last line in the file. To start on a
% new line instead, the input text should contain a newline character at
% the beginning (or a cell with an empty char, if it is a cell array).
%
%
% Written by Wilfried Beslin
% Last Updated 2023-12-01 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    p = inputParser();
    
    p.addRequired('txt', @(v) ischar(v) || iscellstr(v))
    p.addRequired('filePath', @(v) validateattributes(v,{'char'},{'row'}))
    p.addOptional('existingFileOpt', 'prompt', @ischar);
    
    p.parse(txt, filePath, varargin{:})
    existingFileOpt = validatestring(lower(p.Results.existingFileOpt), {'prompt','overwrite','append','cancel'});
    
    % convert cellstr to char
    if iscellstr(txt)
        txt = strjoin(txt, '\n');
    end
    
    % handle existing files
    if isfile(filePath)
        processingExistingFile = true;
        while processingExistingFile
            switch existingFileOpt
                case 'prompt'
                    questPrompt = sprintf('The file "%s" already exists.\n\nWhat would you like to do?', filePath);
                    promptOptOverwrite = 'Overwrite File';
                    promptOptAppend = 'Append to File';
                    promptOptCancel = 'Cancel';
                    existingFileOptLong = questdlg(questPrompt, 'File Exists', promptOptOverwrite, promptOptAppend, promptOptCancel, promptOptCancel);
                    
                    switch existingFileOptLong
                        case promptOptOverwrite
                            existingFileOpt = 'overwrite';
                        case promptOptAppend
                            existingFileOpt = 'append';
                        case promptOptCancel
                            existingFileOpt = 'cancel';
                        otherwise
                            existingFileOpt = '';
                    end
                    
                case 'overwrite'
                    openType = 'w';
                    processingExistingFile = false;
                    
                case 'append'
                    openType = 'a';
                    processingExistingFile = false;
                    
                case 'cancel'
                    errMsg = 'File already exists; operation cancelled';
                    disp(errMsg)
                    return
                    
                otherwise
                    errMsg = 'Operation cancelled';
                    disp(errMsg)
                    return
            end
        end
    else
        openType = 'w';
    end
    
    % try writing file
    try
        [fid, errMsg] = fopen(filePath, openType);
        if fid ~= -1
            fprintf(fid, txt);
            fclose(fid);
        end
    catch ME
        errMsg = ME.message;
    end
    
    % display warning if any errors present
    if ~isempty(errMsg)
        warning(errMsg)
    end
    
end