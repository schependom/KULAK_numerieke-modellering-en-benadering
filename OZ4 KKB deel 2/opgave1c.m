%%% Exercise 1cb

%% Compute integral for different degrees

% Show enough digits precision
format long 

% Number of evaluation points
N = 1000; 

% x-values
x = linspace(1/2,3/2,N)'; 

% Function values
f = expint(1) - expint(x); 

Mmax = 3;

% Compute integral
for M = 1:Mmax
    disp(M);
    disp(integraal(x, f, M));
end