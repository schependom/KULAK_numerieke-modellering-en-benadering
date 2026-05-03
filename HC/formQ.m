function Q = formQ(W)
% FORMQ
%   Computes the matrix Q from the full QR factorization A = QR,
%   given the implicit representation W produced by house(A).
%
%
%       house(A) stores the normalised reflector vectors v_k in the
%       columns of W.  The full orthogonal factor is:
%
%           Q* =  H_n * ... * H_2 * H_1       (product of n reflectors)
%
%       where each H_k = I - (2 * v_k * v_k')
%
%       The reason we don't form H_k has to do with memory and calculation
%       complexity. Instead, we only store the reflectors v_k in the
%       subdiagonal matrix W.
%
%
%       To form Q = QI explicitly, we can compute Qx via Algorithm 10.3,
%       with x being unit vectors e_j in the identity matrix I.
%
%
%   INPUTS
%
%       W   (m x n) lower-triangular matrix from house(A)
%           k-th column holds v_k in positions k:m
%
%   OUTPUTS
%
%       Q   (m x m) orthogonal matrix
%           satisfies A = Q * R  when combined with house's R

    [m, n] = size(W);

    % Start from the (m x m) identity and accumulate reflectors
    Q = eye(m);

    % Apply reflectors in reverse order: H_n, H_{n-1}, ..., H_1
    %
    %   At step k we apply H_k to the trailing block Q(k:m, k:m):
    %
    %       H_k * B  =  B - 2 * v_k * (v_k' * B)
    %
    %   The leading k-1 rows and columns of Q are untouched by H_k,
    %   because v_k has zeros in positions 1:k-1 (W is lower-triangular).

    % ALGORITHM 10.3

    % Begin VANBINNEN!
    %
    % Qx = (Q_1 * ... * Q_k)x
    %    = (Q_1 (Q_2 (... (Q_k x))))
    %
    % Vandaar de k DOWN TO 1, want we BEGINNEN VANBINNEN

    % For k=n down to 1
    for k = n:-1:1

        % v_k = W(n:m, n) = laatste kolom vanaf diagonaal tot beneden
        v = W(k:m, k);

        % Skip zero columns (corresponds to a skipped step in house)
        if norm(v) < eps
            continue;
        end
        
        % x_k:m = x_k:m - 2 v_k (v_k^* x_k:m)
        % e_k:m = e_k:m - 2 v_k (v_k^* e_k:m)

        % Apply H_k = I - 2*v*v' to the trailing block
        Q(k:m, k:m) = Q(k:m, k:m) - 2 * v * (v' * Q(k:m, k:m));
    end

end