function T = tridiag(A)
% TRIDIAG
%
%   Reduces a real symmetric matrix to tridiagonal form
%   by orthogonal similarity transformations.
%
%       This is STEP 1 in the two-step unitary diagonalization 
%       of a real symmetric matrix:
%
%       A=A^T   --STEP 1->  (Q^T)    T    Q
%               --STEP 2->  (Q^T) \Lambda Q
%
%       The Q's in this first step are Householder reflections that
%       keep one dimension untouched
%           -> Q leaves the first column untouched
%           -> Q^T leaves the first row untouched
%       
%       These HH reflections introduces ZEROS beneath the SECOND DIAGONAL
%       (>< classical HH is zero beneath the FIRST diagonal)
%
%           => HESSENBERG FORM
%           => TRIDIAGONAL T in case of A=A^T \in \R^(m x m)
%
%
%   We only use ELEMENTARY Matlab operations (e.g. NO hess())
%
%
%   INPUTS
%       A   Real symmetric matrix
%
%   OUTPUTS
%       T   Tridiagonalization of A, symmetric and tridiagonal
%           up to rounding errors.
    
    % Check whether A is real
    if ~isreal(A)
        error('tridiag: input matrix A must be real.');
    end
 
    % Check whether A is symmetric
    if ~issymmetric(A)
        error('tridiag: input matrix A must be symmetric.');
    end

    % Dimension
    m = size(A,1);

    % m-2 Householder steps to obtain tridiagonal form
    for k = 1:m-2
 
        % Extract the sub-column to be zeroed out.
        % 
        %       x is the part of column k that lies 
        %       strictly below the SUBdiagonal.
        % 
        %       After the reflection we want: x -> ||x|| * e1
        %
        x = A(k+1:m, k);          % (m-k) x 1 vector
 
        % Unit vector e1 of matching size
        e1 = zeros(m-k, 1);
        e1(1) = 1;
 
        % Reflector vector v
        % 
        %       The sign choice avoids catastrophic 
        %       cancellation when x(1) > 0.
        % 
        %       sign(0) returns 0 in MATLAB, so we guard with max(...,1).
        %
        sgn = sign(x(1));
        if sgn == 0
            sgn = 1;
        end

        % v = x + ||x|| * e1
        v = x + sgn * norm(x) * e1;  
 
        % Skip if x is already zero (nothing to reflect)
        if norm(v) < eps
            continue;
        end
 
        % Normalise so that v'*v = 1  =>  H = I - 2*v*v'
        v = v / norm(v);
 
        % Put the (m-k) x (m-k) sub-reflector into an m x m identity
        Qstar = eye(m);
        Qstar(k+1:m, k+1:m) = eye(m-k) - 2 * (v * v'); % H = I - 2*v*v'
 
        % Apply the orthogonal similarity transformation
        % A <- (Q^*) * A * Q
        A = Qstar * A * Qstar';
    end
 
    T = A;

end