t = [1.01; 1.02; 1.03; 1.04; 1.05; 1.06; 1.07];

% monomiale basis
A1 = [ones(7,1) t]
% -> kolommen zijn hoog-collineair
% -> wijzen in zelfde richting

% zorgvuldig gekozen basis
A2 = [ones(7,1) 30*(t-1.04)]
% -> kolommen loodrecht op elkaar
% -> a_21 \perp a_22
% -> goede conditie

cond1 = cond(A1) % grote  kappa = 104.1
cond2 = cond(A2) % kleine kappa = 1.667




