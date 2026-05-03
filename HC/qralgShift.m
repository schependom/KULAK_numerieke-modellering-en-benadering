function [Tnew, residuals] = qralgShift(T)
% QRALGSHIFT
%   Runs the SHIFTED QR ALGORITHM on a real symmetric tridiagonal matrix T,
%   using the WILKINSON SHIFT at each step.
%
%
%   WILKINSON SHIFT
%
%       At each step we shift by mu, chosen as the eigenvalue of the
%       bottom-right 2x2 block of T that is CLOSEST to T(m,m):
%
%           B = [ T(m-1,m-1)   T(m-1,m) ]
%               [ T(m,m-1)     T(m,m)   ]
%
%       The two eigenvalues of B are:
%
%           delta = (T(m-1,m-1) - T(m,m)) / 2
%           mu    = T(m,m) - sign(delta)*T(m,m-1)^2 / (|delta| + sqrt(delta^2 + T(m,m-1)^2))
%
%       This shift makes the (m, m-1) entry converge CUBICALLY (generically).
%
%
%   ALGORITHM
%
%       T_shifted = T - mu*I
%       [Q, R]    = qr(T_shifted)
%       T_new     = R*Q + mu*I          % shift back
%
%   This is equivalent to one step of:
%       T_new = Q'*T*Q                  % same eigenvalues
%
%
%   STOPPING CRITERION
%
%       |t_{m,m-1}| < tol,   tol = 1e-11
%
%
%   INPUTS
%
%       T   Real symmetric tridiagonal (m x m) matrix
%
%   OUTPUTS
%
%       Tnew        Tridiagonal matrix after convergence of t_{m,m}
%       residuals   Vector of |t_{m,m-1}| values per iteration

    m         = size(T, 1);
    tol       = 1e-11;
    residuals = [];
    iter      = 0;

    while abs(T(m, m-1)) >= tol

        % ---- Wilkinson shift ----
        %
        %   delta = (a - d) / 2  where a = T(m-1,m-1), d = T(m,m)
        %   b     = T(m, m-1)
        %
        a     = T(m-1, m-1);
        d     = T(m,   m  );
        b     = T(m,   m-1);

        delta = (a - d) / 2;

        % Eigenvalue of the 2x2 corner closest to T(m,m)
        if delta == 0
            mu = d - abs(b);
        else
            mu = d - sign(delta) * b^2 / (abs(delta) + sqrt(delta^2 + b^2));
        end

        % ---- Shifted QR step ----
        Tshift = T - mu * eye(m);

        [W, R] = house(Tshift);
        Q      = formQ(W);

        T = R * Q + mu * eye(m);   % shift back

        % Enforce symmetry and tridiagonality
        T = (T + T') / 2;
        d_  = diag(T);
        s_  = diag(T, 1);
        T   = diag(d_) + diag(s_, 1) + diag(s_, -1);

        iter = iter + 1;
        residuals(end+1) = abs(T(m, m-1)); %#ok<AGROW>
    end

    fprintf('qralgShift converged in %d iterations.\n', iter);

    Tnew = T;
end