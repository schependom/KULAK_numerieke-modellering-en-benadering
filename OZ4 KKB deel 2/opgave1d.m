%%% Exercise 1d

%% Compute residue of differential approximation for different degrees

% Set up x and f values
format long % Show enough digits precision
N = 1000; % Number of evaluation points
x = linspace(1/2,3/2,N)'; 

% Function values
f = expint(1) - expint(x); 

% Evalue analytic differential
ftrue = exp(-x)./x; 

Mmax = 15;
reslist = [];

% Loop over degrees
for M = 1:Mmax
   disp(M);

   % Calculate  differentials in the points x
   a = afgeleide(x, f, M);

   % Norm of residue (relative)
   res = max(abs(ftrue-a)./abs(ftrue)); 

   % Store norm of residue
   reslist=[reslist;res];

   disp(res);
end

disp(reslist);
semilogy(1:Mmax,reslist)

% The derivative error starts relatively large at low M, 
% decreases rapidly until roughly M \approx 8-10, 
% then flattens or even increases slightly for M>12.
%
% The minimum error is \approx 10e-5
%
% Waarom? TODO

%% Compare to integral

intreslist = [];
Itrue = -0.034159141390239;

for M = 1:Mmax
   disp(M);
   I = integraal(x, f, M);

   % Norm of residue
   intres = max(abs(Itrue-I)./abs(Itrue));

   % Save
   intreslist=[intreslist;intres];

   % Display residu
   disp(intres);
end

semilogy(1:Mmax, intreslist)
% De fout daalt monotoon en bereikt de machineprecisie
%
% WAAROM? TODO!

%% Side by side
disp([reslist,intreslist]);