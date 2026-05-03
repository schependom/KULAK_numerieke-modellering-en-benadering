%%%%%%%%%%%%%%%%%%%%%%
% PAGINA 223-224, OEFENING 29.d
%

%% e)
%
%   Matrix: A = diag(15:-1:1) + ones(15,15)
%
%   Compare UNSHIFTED vs SHIFTED QR on this matrix.
%   Generate two sawtooth plots and discuss convergence.
 
A = diag(15:-1:1) + ones(15, 15);
m = size(A, 1);
 
% Symmetry check
assert(norm(A - A', 'fro') < 1e-12, 'A is not symmetric');
 
T = tridiag(A);
 
tol = 1e-11;
 
%% Helper: deflation driver (returns all residuals)
function [eigenvalues, all_residuals, pass_lengths] = deflate(T, use_shift)
    m             = size(T, 1);
    eigenvalues   = zeros(m, 1);
    all_residuals = [];
    pass_lengths  = [];
    Tcurr         = T;
    tol_          = 1e-11;
 
    for sz = m : -1 : 2
        Tsub = Tcurr(1:sz, 1:sz);
 
        if use_shift
            [Tsub, residuals] = qralgShift(Tsub);
        else
            % inline unshifted qralg that returns residuals
            residuals = [];
            iter_     = 0;
            while abs(Tsub(sz, sz-1)) >= tol_
                [W, R] = house(Tsub);
                Q_     = formQ(W);
                Tsub   = R * Q_;
                Tsub   = (Tsub + Tsub') / 2;
                d_     = diag(Tsub);  s_ = diag(Tsub,1);
                Tsub   = diag(d_) + diag(s_,1) + diag(s_,-1);
                iter_  = iter_ + 1;
                residuals(end+1) = abs(Tsub(sz, sz-1)); %#ok<AGROW>
            end
        end
 
        eigenvalues(m - sz + 1) = Tsub(sz, sz);
        all_residuals            = [all_residuals, residuals]; %#ok<AGROW>
        pass_lengths(end+1)      = numel(residuals); %#ok<AGROW>
 
        Tcurr = Tsub(1:sz-1, 1:sz-1);
    end
    eigenvalues(m) = Tcurr(1,1);
end
 
%% Run both
fprintf('=== Unshifted QR ===\n');
[ev_unshifted, res_unshifted, pl_unshifted] = deflate(T, false);
 
fprintf('\n=== Wilkinson-Shifted QR ===\n');
[ev_shifted,   res_shifted,   pl_shifted  ] = deflate(T, true);
 
lambda_eig = sort(eig(A), 'descend');
 
fprintf('\nMax error unshifted: %.4e\n', max(abs(sort(ev_unshifted,'descend') - lambda_eig)));
fprintf('Max error shifted:   %.4e\n', max(abs(sort(ev_shifted,  'descend') - lambda_eig)));
 
fprintf('\nTotal QR iterations unshifted: %d\n', sum(pl_unshifted));
fprintf('Total QR iterations shifted:   %d\n', sum(pl_shifted));
 
%% Plot 1: Unshifted
figure('Name', 'e) Unshifted QR – diag(15:-1:1)+ones(15)');
semilogy(res_unshifted, 'b-o', 'MarkerSize', 2, 'LineWidth', 1);
hold on;
yline(tol, 'r--', 'Label', 'tol');
cumlen = cumsum(pl_unshifted);
for i = 1:length(cumlen)-1
    xline(cumlen(i)+0.5, 'k:', 'LineWidth', 0.8);
end
xlabel('Cumulative QR iteration');  ylabel('|t_{m,m-1}|');
title('Oef. 29.1e – UNSHIFTED QR on diag(15:-1:1)+ones(15)');
grid on;
saveas(gcf, 'sawtooth_e_unshifted.png');
 
%% Plot 2: Shifted
figure('Name', 'e) Shifted QR – diag(15:-1:1)+ones(15)');
semilogy(res_shifted, 'r-o', 'MarkerSize', 2, 'LineWidth', 1);
hold on;
yline(tol, 'b--', 'Label', 'tol');
cumlen = cumsum(pl_shifted);
for i = 1:length(cumlen)-1
    xline(cumlen(i)+0.5, 'k:', 'LineWidth', 0.8);
end
xlabel('Cumulative QR iteration');  ylabel('|t_{m,m-1}|');
title('Oef. 29.1e – WILKINSON-SHIFTED QR on diag(15:-1:1)+ones(15)');
grid on;
saveas(gcf, 'sawtooth_e_shifted.png');
 
fprintf('Plots saved.\n');
 