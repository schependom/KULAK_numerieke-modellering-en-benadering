%% Givens Transformatie

% vector
x = randn(2,1);

% r = \sqrt{a^2+b^2}
r = norm(x);

% COS = A/S
c = x(1)/r;
% SIN = O/S
s = -x(2)/r; % negatief!

% STANDAARD rotatie TEGEN DE KLOK IN met hoek alpha
% gebeurt als volgt:
%
%       R = [   c           -sin(alpha); 
%               sin(alpha)  c           ]

% Om de b-component op nul te krijgen, moeten we de vector
% echter MET DE KLOK MEE roteren (dus -alpha)!
%
%   ->  Voor een rotatie met -alpha gebruiken 
%       we de getransponeerde matrix G'.
%
%   ->  Daarom dat we bij Givens steeds linksvermenigvuldigen
%       met G^T
%
%   ->  Dit is equivalent aan S=-sin stellen, omdat
%       -sin(alpha) = sin(-alpha).

% Rotatie over alpha MET DE KLOK MEE
G = [c -s; s c];

% (a,b) is gedraaid naar (r,0)
G*x