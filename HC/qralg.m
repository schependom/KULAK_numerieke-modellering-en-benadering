function Tnew = qralg(T)
% QRALG
%   Runs the UNSHIFTED QR ALGORITHM on a real symmetric tridiagonal
%   matrix T until convergence of the bottom eigenvalue.
%
%       The QR algorithm iterates the similarity transformation:
%
%           T_0         =  T
%           Q_k * R_k   =  T_{k-1}
%           T_k         =  R_k * Q_k     (recombined in reverse order)
%
%
%       Herinner: bij Simultaneous Iteration gebruiken we Q_underscore,
%       die hier expliciet zou kunnen uitgerekend worden als het lopende
%       product van de Q's:
%
%           Q_underscore(k) = Q(1) * Q(2) * ... * Q(k)
%
%       Als we nu schrijven
%
%           R_underscore(k) = R(k) * ... * R(2) * R(1),
%
%       Dan krijgen we dat 
% 
%           Q_underscore(k) * R_underscore(k) = A(k)
%
%
%       Aangezien Q_underscore bij simultane iteratie convergeert
%       naar een matrix met in de kolommen de eigenvectoren,
%       kunnen we concluderen dat de similariteitstransformaties
%
%           Q_underscore(k)^T A Q(k)_underscore
%           = [Q(1) * ... * Q(k)]^T Q [Q(1) * ... * Q(k)]
%
%       zullen convergeren naar een diagonaalmatrix met de eigenwaarden
%       op de hoofddiagonaal (releigh coefficienten!)
%
%
%       Convergence is monitored via the (m, m-1) sub-diagonal entry:
%       when it is small enough, the last diagonal entry t_{m,m} has
%       converged to an eigenvalue and we stop.
%       This last diagonal entry will converge the latest.
%
%
%   STOPPING CRITERION
%
%       |t_{m,m-1}| < tol,   tol = 1e-11   (= 10e-12 as in the exercise)
%
%
%   SYMMETRY AND TRIDIAGONALITY ENFORCEMENT
%
%       In exact arithmetic each T_k remains symmetric and tridiagonal.
%       Floating-point errors accumulate over iterations, so after each
%       step we:
%           (1) symmetrize:      T <- (T + T') / 2
%           (2) tridiagonalize:  zero out all entries beyond the first
%                                super- and sub-diagonal
%
%
%   INPUTS
%
%       T   Real symmetric tridiagonal (m x m) matrix
%
%   OUTPUTS
%
%       Tnew    Tridiagonal matrix after convergence of t_{m,m}
%               Diagonal entries approximate the eigenvalues of T.

    m   = size(T, 1);
    tol = 1e-11;

    iter = 0;

    % QR iteration
    %
    %   Each iteration is one QR step:
    %       T = Q*R  =>  T_new = R*Q  (same eigenvalues, more diagonal)
    %
    while abs(T(m, m-1)) >= tol

        % QR factorization via Householder reflections (implicit)
        [W, R] = house(T);

        % Form Q explicitly from the stored reflector vectors
        Q = formQ(W);

        % Recombine: T_new = R * Q
        %
        %   Note: T_new = R*Q = (Q'*T)*Q = Q'*T*Q
        %
        T = R * Q;

        % Enforce symmetry
        %
        %   T should be symmetric after each step (A=A^T, orthog. simil.)
        %   Symmetrize to prevent asymmetric drift under rounding.
        %
        T = (T + T') / 2;

        % Enforce tridiagonality
        %
        %   Extract only the three relevant diagonals and zero the rest.
        %   diag(v, k): k=0 main, k=1 super, k=-1 sub.
        %
        d  = diag(T);           % main diagonal
        s  = diag(T,  1);       % super-diagonal
        T  = diag(d) + diag(s, 1) + diag(s, -1);
        %   sub-diagonal = super-diagonal (symmetry enforced above)

        iter = iter + 1;
    end

    fprintf('qralg converged in %d iterations.\n', iter);

    Tnew = T;

end