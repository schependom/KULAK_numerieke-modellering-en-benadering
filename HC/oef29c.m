%%%%%%%%%%%%%%%%%%%%%%
% PAGINA 223-224, OEFENING 29.c
%

%% c)
% Write a driver program which 
%
%    (i) calls tridiag
%   (ii) calls qralg to get one eigenvalue
%  (iii) calls qralg with a smaller matrix to get another eigenvalue
%   (iv) ... and so on until all of the eigenvalues of A are determined
%
% Set things up so that the values of |t_m,m-1| at every QR iteration
% are stored in a vector and so that at the end, the program generates a
% semilogy plot of these values as a function of the number of QR
% factorizations.
%
% Here m will step from length(A) to length(A)-1 and so on down to 3 and
% finally 2 as the deflation proceeds, and the plot will be correspondingly
% sawtoothed.
%
% Run your program for A = hilb(4). The output should be a set of
% eigenvalues and a sawtooth plot.

A = hilb(4);
m = size(A, 1);

%%
% Tridiagonalisation
T = tridiag(A);

%%
% Deflation loop
%
%   We shrink the active submatrix after each eigenvalue converges:
%
%       Pass 1:  full  4x4  T,  find lambda_4  (smallest)
%       Pass 2:  top-left 3x3 block,  find lambda_3
%       Pass 3:  top-left 2x2 block,  find lambda_2
%       Pass 4:  1x1 block  ->  lambda_1 trivially
%
eigenvalues   = zeros(m, 1);
all_residuals = [];          % will grow: concatenation of each pass
pass_lengths  = [];          % number of iters per pass, for x-axis marks
 
tol = 1e-11;
 
Tcurr = T;
 
for sz = m : -1 : 2
 
    Tsub = Tcurr(1:sz, 1:sz);
 
    % ---- inline qralg that also records residuals ----
    residuals = [];
    iter = 0;
 
    while abs(Tsub(sz, sz-1)) >= tol
 
        [W, R] = house(Tsub);
        Q      = formQ(W);
        Tsub   = R * Q;
 
        % enforce symmetry and tridiagonality
        Tsub = (Tsub + Tsub') / 2;
        d    = diag(Tsub);
        s    = diag(Tsub, 1);
        Tsub = diag(d) + diag(s,1) + diag(s,-1);
 
        iter = iter + 1;
        residuals(end+1) = abs(Tsub(sz, sz-1)); %#ok<AGROW>
    end
 
    fprintf('Pass (sz=%d): converged in %d iterations. lambda_%d = %.10f\n', ...
            sz, iter, m-sz+1, Tsub(sz,sz));
 
    eigenvalues(m - sz + 1) = Tsub(sz, sz);
    all_residuals            = [all_residuals, residuals]; %#ok<AGROW>
    pass_lengths(end+1)      = iter; %#ok<AGROW>
 
    % deflate: keep the top (sz-1)x(sz-1) block for the next pass
    Tcurr = Tsub(1:sz-1, 1:sz-1);
end
 
% Last eigenvalue is the remaining 1x1
eigenvalues(m) = Tcurr(1,1);
fprintf('Pass (sz=1): trivial. lambda_%d = %.10f\n', m, eigenvalues(m));
 
%%
% Sort and compare
lambda_qralg = sort(eigenvalues, 'descend');
lambda_eig   = sort(eig(A),      'descend');
 
fprintf('\nEigenvalues (qralg):  '); fprintf('%.8f  ', lambda_qralg); fprintf('\n');
fprintf('Eigenvalues (eig):    '); fprintf('%.8f  ', lambda_eig);   fprintf('\n');
fprintf('Max error: %.4e\n', max(abs(lambda_qralg - lambda_eig)));

%%
% Sawtooth semilogy plot
figure('Name', 'c) Unshifted QR – sawtooth');
semilogy(all_residuals, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
hold on;
yline(tol, 'r--', 'LineWidth', 1, 'Label', 'tol = 1e-11');
 
% Mark pass boundaries with vertical lines
cumlen = cumsum(pass_lengths);
for i = 1:length(cumlen)-1
    xline(cumlen(i) + 0.5, 'k:', 'LineWidth', 1);
end
 
xlabel('QR iteration (cumulative)');
ylabel('|t_{m,m-1}|');
title('Oef. 29.1c – Unshifted QR on hilb(4): sawtooth convergence');
grid on;
legend('residual |t_{m,m-1}|', 'tolerance', 'Location', 'southwest');
 
saveas(gcf, 'sawtooth_unshifted.png');
fprintf('Plot saved.\n');