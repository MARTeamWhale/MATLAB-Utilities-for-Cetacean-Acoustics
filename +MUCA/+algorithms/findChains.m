function [chainStart, chainStop, chainLength] = findChains(vec, target, varargin)
%
% Look for successive repetitions of a certain value in a vector.
%
% SYNTAX:
%   [chainStart, chainStop, chainLength] = findChains(vec, target)
%   [__] = findChains(vec, target, nMin)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "vec" - vector of values within which to look for repetitions
%   .......................................................................
%   "target" - the value to look for
%   .......................................................................
%   "nMin" [OPTIONAL] - minimum chain length to detect. Default is 1.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "chainStart" - start index of each chain
%   .......................................................................
%   "chainStop" - stop index of each chain
%   .......................................................................
%   "chainLength" - number of elements in each chain
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
% Last Updated 2022-03-22 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    p = inputParser();
    p.addRequired('vec', @isvector)
    p.addRequired('target', @isscalar)
    p.addOptional('nMin', 1, @(v) validateattributes(v, {'numeric'},{'scalar','integer','positive'}))
    
    p.parse(vec, target, varargin{:})
    nMin = p.Results.nMin;

    % determine if row or column vector
    if isrow(vec)
        dim = 2;
    else
        dim = 1;
    end
    n = numel(vec);
    
    % find hits
    hit = vec == target;
    
    % convert hits to +1/-1
    hitPolar = int8(hit);
    hitPolar(hitPolar == 0) = -1;
    
    % find zero crossings in polar hit vector
    zc = false(size(hitPolar));
    zc(2:end) = hitPolar(1:end-1).*hitPolar(2:end) < 0; % Definition of zero crossing events
    
    % determine when events start and end
    chainStart = find(zc & hit);
    chainStop = find(zc & ~hit) - 1;
    
    if chainStop(1) < chainStart(1)
        chainStart = cat(dim, 1, chainStart);
    end
    if chainStart(end) > chainStop(end)
        chainStop = cat(dim, chainStop, n);
    end
    
    % determine chain length
    chainLength = chainStop - chainStart + 1;
    
    % trim results
    keep = chainLength >= nMin;
    chainStart = chainStart(keep);
    chainStop = chainStop(keep);
    chainLength = chainLength(keep);
end