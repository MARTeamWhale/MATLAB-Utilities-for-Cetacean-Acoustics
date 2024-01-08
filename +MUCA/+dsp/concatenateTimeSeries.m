function [xCat,tCat] = concatenateTimeSeries(xPieces, fs, tStarts, varargin)
%
% Pieces together multiple timeseries vectors into a single vector, keeping
% the sampling rate consistent and conserving the timestamps of the
% original vectors as much as possible. This means the final vector could
% have samples shifted in time very slightly to ensure a consistent time
% delta. If there are gaps between files, they will be filled in with
% either zeros or NaNs (depending on user preference). If there are
% timeseries that overlap, then some will have their first couple of
% samples removed - a warning will be issued in this case.
%
% Currently only works with single-channel data.
%
% SYNTAX:
%   [xCat, tCat] = concatenateTimeSeries(xPieces, fs, tStarts)
%   [xCat, tCat] = concatenateTimeSeries(xPieces, fs, tStarts, 'PadType',padTypeStr)
%
% INPUT ARGUMENTS:
%   Required
%   .......................................................................
%   "xPieces" - Cell array containing 1-D time series vectors
%   .......................................................................
%   "fs" - Sampling rate
%   .......................................................................
%   "tStarts" - Vector specifying the start time of each time series piece.
%       May be numeric, Duration, or Datetime. Must be monotonic.
%   .......................................................................
%
%   Optional (Name-Value pairs)
%   .......................................................................
%   "PadType" - Char string specifying the type of data to use to fill in
%       gaps. Options are:
%           'zeros' (default)
%           'nans'
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "xCat" - Concatenated time series samples
%   .......................................................................
%   "tCat" - Times of concatenated time series samples
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
% Last updated 2024-01-08 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    p = inputParser();
    
    % input parsing and validation
    p.addRequired('xPieces', @iscell)
    p.addRequired('fs', @(val)validateattributes(val,{'numeric'},{'scalar','positive'}));
    p.addRequired('tStarts', @(val)validateattributes(val,{'numeric','duration','datetime'},{'vector'}));
    p.addParameter('PadType','zeros',@(val)validateattributes(val,{'char'},{'vector'}));
    
    p.parse(xPieces, fs, tStarts, varargin{:});
    fs = double(fs);
    padType = validatestring(p.Results.PadType, {'zeros','nans'});
    
    numPieces = numel(xPieces);
    assert(numel(tStarts) == numPieces, 'The number of elements in xPieces and tStarts must match!');
    
    % initialize variables
    delta_t_fs = 1/fs;
    delta_t_tol = delta_t_fs/2;
    tAbsStart = tStarts(1);
    tRelStarts = tStarts - tAbsStart;
    usingTimeObj = isduration(tRelStarts);
    if usingTimeObj
        tRelStarts = seconds(tRelStarts);
    end
    switch padType
        case 'zeros'
            padFcn = @zeros;
        case 'nans'
            padFcn = @NaN;
        otherwise
            error('Unrecognized value for "PadType"')
    end
    hasOverlap = false;
    
    % process each piece
    for ii = 1:numPieces
        x_ii = xPieces{ii};
        tStart_ii = tRelStarts(ii);
        
        % ensure column vectors
        if isrow(x_ii)
            x_ii = x_ii';
        elseif ~iscolumn(x_ii)
            error('Data pieces must be vectors!')
        end
        
        % create relative time vector for this piece
        numSamples_ii = numel(x_ii);
        tRel_ii = (0:(numSamples_ii-1))'/fs;
        
        if ii == 1
            % process first piece
            
            % get and validate data type
            assert(isnumeric(x_ii), 'Data must be numeric!')
            xType = class(x_ii);
            
            % initialize output
            xCat = x_ii;
            tRelCat = tRel_ii;
            
        else
            % process remaining pieces
            
            % ensure the type is correct
            assert(strcmp(xType,class(x_ii)), 'Data types are not consistent!')
            
            % get amount of separation
            delta_t_sep = tStart_ii - tRelCat(end);
            
            % handle each separation case - this is where the magic happens
            if abs(delta_t_sep-delta_t_fs) <= delta_t_tol
                % Separation is within tolerance - shift the vector being
                % added as needed
                xAppend_ii = x_ii;
                tAppend_ii = tRelCat(end) + delta_t_fs + tRel_ii;
                
            elseif delta_t_sep - delta_t_fs > delta_t_tol
                % Separation is positively greater than tolerance, creating
                % a gap; in this case, insert padding as necessary
                nZeros_ii = fix(delta_t_sep/delta_t_fs);

                xAppend_ii = [feval(padFcn,nZeros_ii,1,xType); x_ii];
                tAppend_ii = tRelCat(end) + [delta_t_fs*(1:nZeros_ii)'; tRel_ii + delta_t_fs*(nZeros_ii+1)];
                
            else
                % Separation is negatively greater than tolerence, creating
                % an overlap; in this case, remove the first extra samples
                % from the vector being added. Issue a warning in this
                % case, because it results in data destruction.
                hasOverlap = true;
                
                tRelShifted_ii = tRel_ii - abs(delta_t_sep);
                nOverlap_ii = sum(tRelShifted_ii < 0);
                
                xAppend_ii = x_ii((nOverlap_ii+1):end);
                tAppend_ii = tRelCat(end) + delta_t_fs + tRel_ii(1:(end-nOverlap_ii));
            end
            
            % concatenate vectors
            xCat = [xCat; xAppend_ii];
            tRelCat = [tRelCat; tAppend_ii];
        end
    end
    
    % create absolute time vector by adding absolute time
    if usingTimeObj
        tRelCat = seconds(tRelCat);
    end
    tCat = tAbsStart + tRelCat;

    % issue warning if overlap was detected
    if hasOverlap
        warning('Some timeseries pieces were overlapping and had to be truncated. Results could contain missing and/or inaccurate data.')
    end
end