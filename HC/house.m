function [W, R] = house(A)
% HOUSE
%   Computes an implicit representation of a full QR factorization
%   A = QR of an (m x n) matrix with (m >= n) using HH reflections.
%
%       At each step k, a Householder reflector H_k is chosen such that
%
%           H_k * A(k:m, k)  =  -sgn(x_1) * ||x|| * e_1
%
%       and applied to the trailing submatrix:
%
%           A(k:m, k:n)  <-  H_k * A(k:m, k:n)
%
%       The reflector vectors v_k are stored in the lower-triangular
%       part of W (columns of W), so Q never needs to be formed
%       explicitly. This is the IMPLICIT representation of Q.
%
%
%   INPUTS
%
%       A   (m x n) matrix with (m >= n)
%
%   OUTPUTS
%
%       W   (m x n) lower-triangular matrix
%           whose k-th column stores the normalised reflector vector v_k
%           in positions k:m  (entries 1:k-1 are zero by construction)
%
%       R   (n x n) upper-triangular matrix
%           upper-triangular part of the reduced A

    [m, n] = size(A);

    % W stores the normalised reflector vectors column by column
    W = zeros(m, n);

    % n Householder steps to reduce A to upper-triangular form
    for k = 1:n

        % Sub-column to be zeroed out below the diagonal
        x = A(k:m, k);          % (m-k+1) x 1

        % Unit vector e1 of matching size
        e1    = zeros(m-k+1, 1);
        e1(1) = 1;

        % Reflector vector v
        %
        %   Sign choice avoids catastrophic cancellation:
        %   we push x AWAY from e1, not toward it.
        %
        %   guard: sign(0) = 0 in MATLAB, default to +1
        %
        sgn = sign(x(1));
        if sgn == 0
            sgn = 1;
        end

        v = x + sgn * norm(x) * e1;

        % Skip if column is already zero (nothing to reflect)
        if norm(v) < eps
            continue;
        end

        % Normalise: ensures H = I - 2*v*v'  with  v'*v = 1
        v = v / norm(v);

        % Store normalised reflector in column k of W
        %
        %   W is lower-triangular: entries 1:k-1 remain zero.
        %
        W(k:m, k) = v;

        % Apply H_k to the trailing submatrix B
        %
        %   H_k * B  =  (I - 2*v*v') * B  =  B - 2*v*(v'*B)
        %
        %   Computing v*(v'*B) exploits the rank-1 structure:
        %   O(m*n) instead of forming H_k explicitly: O(m^2).
        %
        A(k:m, k:n) = A(k:m, k:n) - 2 * v * (v' * A(k:m, k:n));
    end

    % Upper-triangular part of the reduced A is R
    % (entries below diagonal are numerically zero up to rounding)
    R = triu(A(1:n, 1:n));

end