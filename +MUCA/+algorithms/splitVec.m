function [chunkRanges, nChunks] = splitVec(n, chunkSize)
%
% Function for breaking a vector into smaller chunks of fixed size. Returns
% indicies that mark the start and end of chunks to break the vector into.
% This is useful, for example, to split up a long audio time series into
% shorter segments.
%
% SYNTAX:
%   [chunkRanges, nChunks] = splitVec(n, chunkSize)
%
% INPUT ARGUMENTS:
%   Required
%   .......................................................................
%   "n" - size of vector to be chunked
%   .......................................................................
%   "chunkSize" - standard size of chunks. If the vector does not divide
%       evenly by this number, then the last chunk will be shorter.
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "chunkRanges" - N-by-2 matrix of chunk interval indices, where N is the
%       number of chunks.
%   .......................................................................
%   "nChunks" - number of chunks
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
% Last updated 2023-11-30 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    iStarts = (1:chunkSize:n)';
    iEnds = [chunkSize:chunkSize:n, n]';
    chunkRanges = [iStarts, iEnds];
    nChunks = size(chunkRanges, 1);
end