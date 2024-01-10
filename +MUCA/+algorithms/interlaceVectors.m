function c = interlaceVectors(a, b)
%
% Given two vectors A and B, where:
%   A = [A1, A2, A3, ... An]
% and
%   B = [B1, B2, B3, ... Bn]
% return a new vector C of length 2n that consists of:
%   C = [A1, B1, A2, B2, A3, B3, ... An, Bn]
%
% A and B may be either row or column vectors, but they must have the same
% dimensions and data type.
%
% SYNTAX:
%   c = interlaceVectors(a, b)
%
% INPUT ARGUMENTS:
%   .......................................................................
%   "a" - first vector
%   .......................................................................
%   "b" - second vector
%   .......................................................................
%
% OUTPUT ARGUMENTS:
%   .......................................................................
%   "c" - interlaced vector that is a combination of "a" and "b"
%   .......................................................................
%
% OUTPUT FILES:
%   <non>
%
% DEPENDENCIES:
%   <none>
%
%
% Written by Wilfried Beslin
% Last updated 2024-01-10, using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    assert(isvector(a), 'Input must be vectors!')
    assert(all(size(a) == size(b)), 'Vectors must be the same size!');

    n = length(a);
    m = n*2;

    c = repelem(a, 2);

    c(2:2:m) = b;
end