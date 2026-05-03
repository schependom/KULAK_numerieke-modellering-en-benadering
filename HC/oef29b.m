%%%%%%%%%%%%%%%%%%%%%%
% PAGINA 223-224, OEFENING 29.1.b
%

%% From a)

A = hilb(4)
T = tridiag(A)

%% b)
% Write a function Tnew = qralg(T) that runs the UNSHIFTED QR ALGORITHM
% on a real tridiagonal matrix T. 
% 
% For the QR factorization at each step, use programs [W, R] = house(A) 
% and Q = formQ(W) of Exercise 10.2 if available, or MATLABs command qr, 
% or for greater efficiency, a new code based on Givens rotations 
% or (2 x 2) Householder reflections.
%
% Again you may wish to enforce symmetry and tridiagonality at each step.
%
% Your program should stop and return the current tridiagonal matrix T as
% Tnew when the (m, m-1) element satisfies |t_{m,m-1}| < 10e-12.
%
% Again, apply your program to A = hilb(4)

% Start from the tridiagonalization computed in part a)
Tnew = qralg(T)
%
% The QR algorithm converges to a diagonal matrix whose entries
% are the eigenvalues of A (= eigenvalues of T, preserved by
% orthogonal similarity).
%
% Expected: diagonal entries ~= eigenvalues of hilb(4).
%
% Verify against MATLAB's built-in eig():
lambda_qralg = sort(diag(Tnew), 'descend')
lambda_eig   = sort(eig(A),     'descend')

% 2.3823e-10
max_error = max(abs(lambda_qralg - lambda_eig))
%
% trace is preserved throughout (eigenvalues sum = trace of A)
%
% trace(A)    = 1.6762
% sum(lambda) = 1.6762
sum(lambda_qralg)
trace(A)



%% d)
% Create a modified function qralgShift that additionally uses the
% Wilkinson shift at each step. 
% 
% Turn in the new sawooth plot for the same example.

A = hilb(4);
m = size(A, 1);
 
T     = tridiag(A);
Tcurr = T;
 
eigenvalues   = zeros(m, 1);
all_residuals = [];
pass_lengths  = [];
 
tol = 1e-11;
 
for sz = m : -1 : 2
 
    Tsub = Tcurr(1:sz, 1:sz);
 
    [Tsub, residuals] = qralgShift(Tsub);
 
    fprintf('Pass (sz=%d): converged in %d iters. lambda_%d = %.10f\n', ...
            sz, numel(residuals), m-sz+1, Tsub(sz,sz));
 
    eigenvalues(m - sz + 1) = Tsub(sz, sz);
    all_residuals            = [all_residuals, residuals]; 
    pass_lengths(end+1)      = numel(residuals);
 
    Tcurr = Tsub(1:sz-1, 1:sz-1);
end
 
eigenvalues(m) = Tcurr(1,1);
 
lambda_qralg = sort(eigenvalues, 'descend');
lambda_eig   = sort(eig(A),      'descend');
 
fprintf('\nEigenvalues (shifted qralg):  '); fprintf('%.8f  ', lambda_qralg); fprintf('\n');
fprintf('Eigenvalues (eig):            '); fprintf('%.8f  ', lambda_eig);   fprintf('\n');
fprintf('Max error: %.4e\n', max(abs(lambda_qralg - lambda_eig)));
 
% Sawtooth plot
figure('Name', 'd) Wilkinson-shifted QR – sawtooth');
semilogy(all_residuals, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2);
hold on;
yline(tol, 'b--', 'LineWidth', 1, 'Label', 'tol = 1e-11');
 
cumlen = cumsum(pass_lengths);
for i = 1:length(cumlen)-1
    xline(cumlen(i) + 0.5, 'k:', 'LineWidth', 1);
end
 
xlabel('QR iteration (cumulative)');
ylabel('|t_{m,m-1}|');
title('Oef. 29.1d – Wilkinson-shifted QR on hilb(4): sawtooth convergence');
grid on;
legend('residual |t_{m,m-1}|', 'tolerance', 'Location', 'southwest');
 
saveas(gcf, '/mnt/user-data/outputs/sawtooth_shifted.png');
fprintf('Plot saved.\n');

%% e)
% Rerun your program for the matrix A=diag(15:-1:1)+ones(15,15) and
% generate two sawtooth plots corresponding to the shift and no shift.
%
% Discuss the rates of convergence observed here and for the earlier
% matrix. Is the convergence lineear, superlineair, quadratic, cubic, ...?
%
% Is it meaningful to speak of a certain "number of QR iterations per EV"?