function varargout = getRecTimes(audioFilePaths)
%
% Get start and end times of audio files.
%
% SYNTAX:
%   dt = getRecTimes(audioFilePaths)
%   [dtStart, dtEnd] = getRecTimes(audioFilePaths)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "audioFilePaths" - cell array of char strings specifying paths to a
%       series of audio files
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "dt" - n-by-2 matrix of datetime objects representing the start and end
%       times of each audio file
%   .......................................................................
%   "dtStart" - vector of datetime objects representing the start time of
%       each audio file
%   .......................................................................
%   "dtSEnd" - vector of datetime objects representing the end time of
%       each audio file
%   .......................................................................
%
% OUTPUT FILES:
%   <none>
%
% DEPENDENCIES:
%   MUCA.time.readDateTime
%   MUCA.time.isoFormat
%
%
% Written by Wilfried Beslin
% Last Updated 2023-11-30 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    import MUCA.time.readDateTime
    import MUCA.time.isoFormat
    
    nargoutchk(1,2);

    numRecs = numel(audioFilePaths);

    dtStart = readDateTime(audioFilePaths, isoFormat('simplified','milliseconds'));
    
    dtStop = NaT(numRecs, 1);
    for ii = 1:numRecs
        filePath_ii = audioFilePaths{ii};
        try
            info_ii = audioinfo(filePath_ii);
            dtStop(ii) = dtStart(ii) + seconds(info_ii.Duration);
        catch ME
            warning('Could not read end time for file "%s":\n%s', filePath_ii, ME.message)
        end
    end
    dtStop.Format = dtStart.Format;
    
    if nargout == 1
        varargout = {[dtStart, dtStop]};
    elseif nargout == 2
        varargout = {dtStart, dtStop};
    end
end