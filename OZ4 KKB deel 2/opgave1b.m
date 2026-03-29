%%% Exercise 1ca

%% Define x-values, function values and weights

N = 5 % Choose number of approximation points
x = linspace(1/2,3/2,N)'; % x-values
f = expint(1) - expint(x); % Function values
w = ones(size(x)); % Weights

%% Compute least-squares polynomial approximation for different degrees

% Set max degree
% Mmax = N-1; % voor een overgedetermineerd stelsel
Mmax = 2*N; % om te laten zien dat het residu nul wordt

res = zeros(1,Mmax);

% Loop over degrees of approdximation
for M = 1:Mmax

    % Approximate using normal equations (with Hankel matrix)
    % -> compute coefficients in monomial basis
    c = kkb2(x,f,w,M);

    % Evaluate polynomial approximation
    r = polyval(c(end:-1:1),x); 

    % Compute residue
    res(M) = norm(f-r); 

end

%% Plot residues for different degrees

close all;
figure;
semilogy(1:Mmax,res);
grid on;