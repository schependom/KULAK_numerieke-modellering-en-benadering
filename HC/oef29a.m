%%%%%%%%%%%%%%%%%%%%%%
% PAGINA 223-224, OEFENING 29.1.a
%

%% a)
% Write a function T = tridiag(A) that reduces a real symmetric (m x m)
% matrix to tridiagonal form by orthogonal similarity transformations.

% H = hilb(N) is the N-by-N matrix with elements 1/(i+j-1), 
% which is a famous example of a badly conditioned matrix.
A = hilb(4)

% Moderately ill-conditioned: 
%   
%   (sigma_1) / (sigma_m) = 15 000
%   We could lose ~4 digits of precision if we solve Ax = b.
%
%   1.5514e+04
cond(A)

T = tridiag(A)
%
% Correct!
%
% [ d1   s1             ]
% [ s1   d2    s2       ]
% [      s2    d3    s3 ]
% [            s3    d4 ]
%
% Symmetry is preserved: the super- and sub-diagonals are equal.
% Off-tridiagonal entries are numerical zero. => tridiagonal

% trace(A) = 1 + 0.3333 + 0.2000 + 0.1429       = 1.6762
% trace(T) = 1.0000 + 0.6506 + 0.0253 + 0.0003  = 1.6762
trace(A)
trace(T)