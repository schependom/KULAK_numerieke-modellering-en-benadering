%%%%%%%%%%%%%%%%%%%%%%
% PAGINA 223-224, OEFENING 29.d
%

%% d)
%
%   Same deflation driver as c), but using the WILKINSON-SHIFTED
%   QR algorithm (qralgShift).
 
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
    all_residuals            = [all_residuals, residuals]; %#ok<AGROW>
    pass_lengths(end+1)      = numel(residuals); %#ok<AGROW>
 
    Tcurr = Tsub(1:sz-1, 1:sz-1);
end
 
eigenvalues(m) = Tcurr(1,1);
 
lambda_qralg = sort(eigenvalues, 'descend');
lambda_eig   = sort(eig(A),      'descend');
 
fprintf('\nEigenvalues (shifted qralg):  '); fprintf('%.8f  ', lambda_qralg); fprintf('\n');
fprintf('Eigenvalues (eig):            '); fprintf('%.8f  ', lambda_eig);   fprintf('\n');
fprintf('Max error: %.4e\n', max(abs(lambda_qralg - lambda_eig)));
 
%% Sawtooth plot
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
 
saveas(gcf, 'sawtooth_shifted.png');
fprintf('Plot saved.\n');