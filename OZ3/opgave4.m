function [A, b, x, kapA, kapbx, eta, costheta] = opgave4(sigma, c, U, V)
%OPGAVE4 Berekent A, b, x en de theoretische waarde voor kappa_b->x
%
% Inputs
%   sigma   vector met n diagonaalelementen
%   c       Rechterlid van S * z = c        (n x 1)
%   U       Orthogonale matrix              (m x m)
%   V       Orthogonale matrix              (n x n)

m = size(U,1); % aantal rijen
n = size(V,1); % aantal rijen

% diagonaalmatrix (n x n)
% en daaronder (m-n) rijen van n nullen
S = [diag(sigma); zeros(m-n,n)]; 

A = U * S * V';
b = U * c;

% \Sigma * z = c <=> z = c ./ (sigma_1, ..., sigma_n)
z = c(1:n) ./ sigma; % elementsgewijze deling

x = V * z;           % vectorverm.

% Conditiegetal van A
kapA = cond(A);

% \eta = ||A||*||x||/||Ax||
eta = norm(A) * norm(x) / norm(A*x);

% cos(theta) = A/S
costheta = norm(A*x) / norm(b);

% Theoretische kappa_b->x
kapbx = kapA / (eta*costheta);
