%%% Kleinste kwadratenproblemen: 
%%% normaalvergelijkingen versus QR-factorisatie

% Merk op dat 
%       A'  = A^H
%       A.' = A^T

% Construeer een voorbeeld waarvoor de normaalvergelijkingen 
% falen en Householder niet.

e = 1e-12;
A = [1 1; 0 e; 0 0];    % is al bovendriehoeks
[Q,R] = qr(A)           % => eenheidsmatrix * A (omdat A = R)
b = rand(3,1);

% Merk op: A^T A = [ 1   1   ]
%                  [ 1 1+e^2 ]
A.' * A
% dit geeft [1 1 ; 1 1], een SINGULIERE matrix zonder inverse!

% De +\epsilon term gaat verloren!!
%
%   1+ (10^{-12})^2 = 1 + 10^{-24}  
%
% => MATLAB maakt hier 1 van!!

%%%
%%% VIA NORMAALVERGELIJKINGEN
%%%

x1 = (A'*A) \ (A'*b) % !!
%
% "x1 is de oplossing van" 
%
%       (A^T A) x = A^T b
%
% MATLAB gebruikt hiervoor:
%       -> Cholesky ontbinding als matrix hermitisch is 
%               M = R^T R
%       -> LU ontbinding als matrix niet hermitisch 
%          of NIET POSITIEF DEFINIET is
%              PM = LU
%
% A^* A is effectief hermitisch (symmetrisch), maar is 
% !! NIET positief definiet !!
% en hier zal dus een LU decompositie gebruikt worden.
%
% Merk op dat (A^T A) ALTIJD PSD (semi!) is
% en bovendien ook strikt PD als de kolommen van A
% onafhankelijk zijn (wat hier door afronding niet meer
% het geval is...)

% Hierboven staan de normaalvergelijkingen:
%
% Residu: r = (b - y_LS)
%           = (b - A x_LS)
%
% Moet loodrecht staat op de kolomruimte van A:
%
%   r \perp C(A) <=> A^T r = 0
%                <=> (A^T A) x = A^T b

% Omdat in (A^T A) de epsilon verloren gaat, wordt deze matrix
% singulier (det=0) en bestaat de inverse dus niet. 
%
% Bij het maken van de LU ontbinding, verschijnt er een nul op de
% diagonaal en kan er geen backsubstitution worden uitvoerd.
%
% Resultaat: x1 = [-Inf ; Inf]

%%%
%%% VIA QR-factorisatie (A = QR)
%%%

x2 = R \ (Q'*b)
%                       A x = b
%                     Q R x = b
% x is oplossing van    R x = Q^T b

% We vermijden hier de berekening van (A^T A) en krijgen een
% hele GROTE resulterende vector: 
%   x2 = [-5.8603e+11 ; 5.8603e+11]
%
% Het is logisch dat x_LS groot is:
%       -> De matrix is BIJNA singulier.
%       -> Kolommen lijken heel hard op elkaar.
%       -> De oplossing moet heel groot zijn om kleine verschil
%          tussen de kolommen te compenseren.

%%%
%%% VERGELIJKING
%%%

% Bereken y_LS - b = Ax_LS - b

norm(A*x1 - b) % NaN
norm(A*x2 - b) % ||r|| = 0.773 = d(b, C(A))