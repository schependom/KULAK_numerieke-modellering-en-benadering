% OEFENING 3
%
% Implementeer het DFT vermenigvuldigingsalgoritme uit de cursus 
% voor twee getallen a en b.
%
% We gebruiken hiervoor de functies num2vec.m
%                                   vec2str.m
%                                   hpmult.m
%
% Met behulp van onze High Precision Multiplication (HP mult) kunnen
% we 100 faculteit exact berekenen. We vergelijken met de normale precisie
% berekening van 100!

%%
% Op pagina 100, 102 wordt beschreven hoe we 2 veeltermen kunnen verm.
% met de DFT + IDFT:
%
%   - we willen p(x) en q(x) vermenigvuldigen met elkaar
%   - zet de graad van het resultaat als N-1
%   - neem als interpolatiepunten x_k = W_n^k
%   - neem DFT van coefficientenrijen van p(x) en q(x)
%   - doe puntsgewijze product hiervan
%   - bereken de IDFT hiervan

%% Compute x!

n = 1;
fac = 100; % Faculty to compute

% Multiply iteratively up to requested number
for k = 2:fac

    % Convert number to character array
    %
    %       e.g. 123 -> [3 2 1]
    vecn = num2vec(k);    

    % Multiply using fft
    n = hpmult(n, vecn);   
    
end

% Convert back to string and print
n = vec2str(n);
fprintf('%i! is equal to %s.\n', fac, n);

%% Gewone precisie
fprintf('%i! in normal precision is equal to %d.\n', fac, gamma(101));