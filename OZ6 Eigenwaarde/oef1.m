%%%%
% A1
%%%%
L1 = diag(1:4); % \Lambda
P1 = orth(rand(4)); % random orthogonale matrix
A1 = P1 * L1 * P1';


%%%%
% A2
%%%%

L2 = L1; % zelfde \Lambda matrix

% Construeer een matrix P2 met conditiegetal K(P2)=100.
%
% Er geldt dat  K(P2) = sigma_1 / sigma_4.
%                 P2  = U S V^T
%
% Hierbij zijn U en V willekeurige orthogonale matrices.
%
% M.a.w. we kiezen S als een diagonaalmatrix waarbij de 
% verhouding tussen de grootste (100) en kleinste (1) waarde 
% precies 100 is (gegenereerd door logspace).
%
P2 = orth(rand(4))*diag(logspace(0, 2, 4))*orth(rand(4))';

% Vermijd expliciet uitrekenen van inverse
% want P2 is NIET orthogonaal, dus P2'!=(P2)^-1
A2 = (P2*L2)/P2; 


%%%%
% A3
%%%%

% Licht aangepaste \Lambda matrix
% zodat A_3 niet meer diagonaliseerbaar is!!
L3 = diag([1 2 3 3]);
L3(3,4) = 1; % This creates a Jordan Block!

% -> we hebben een eigenwaarde herhaald (3)
%
% -> omdat we nu de 1 plaatsen op (3,4), hebben we een Jordan
%    Block gecreeerd, waardoor de geometrische multipliciteit
%    (1) kleiner is dan de algebraische multipliciteit (2)
%
% -> maar er moet gelden dat GM = AM en niet GM < AM!!

% Zelfde P
P3 = P2;

% Vermijd expliciet uitrekenen van inverse
% want P3 is NIET orthogonaal, dus P3'!=(P3)^-1
A3 = (P3*L3)/P3;

%%
% Eerste bepalen we de eigenwaarden met MATLAB.

ew1 = eig(A1)' % 1, 4, 2, 3     -> oke
ew2 = eig(A2)' % 1, 2, 4, 3     -> oke
ew3 = eig(A3)' % 1, 2, 3, 3     -> NOK! Wegens de perturbatie... 

%%
% nlaqr.m berekent met de QR-methode een eigenwaarde van de matrix A.
%
%       e   = gevonden eigenwaarde
%       res = normen voor elke stap

[e1, res1] = nlaqr(A1, true); % true -> pauzes afzetten
                              %      -> toon alle stappen in 1 keer
[e2, res2] = nlaqr(A2, true);
[e3, res3] = nlaqr(A3, true);

% Berekende eigenwaarden
e1
e2
e3

% Plot de convergentie
n1 = length(res1); n2 = length(res2); n3 = length(res3);
semilogy(1:n1, res1, '+-', 1:n2, res2, 'x-', 1:n3, res3, '.-')
xlabel('Iteraties');
ylabel('Residus');
legend('A1','A2','A3');

%%%
% CONVERGENTIE
%
% A3:       -> traagst
%           -> eerder lineaire
%
% A1 & A2:  -> snelst
%           -> moeilijk te zien hoe snel precies
%           -> andere plot nodig

%%%
% BEPALEN VAN DE CONVERGENTIE-ORDE
%
% In de vorige plot (semilogy) zagen we DAT de residus kleiner worden, 
% maar de exacte 'snelheid' (de orde van convergentie) is daar lastig 
% uit af te leiden. 
%
% Om de orde te bepalen, plotten we het huidige residu (res_{k+1}) 
% ten opzichte van het vorige residu (res_k) in een log-log grafiek.
%
% WAAROM KUNNEN WE HIEROP CONVERGENTIE AFLEZEN?
% 
%   De algemene formule voor de fout bij convergentie is:
%       res_{k+1} ≈ C * (res_k)^p
%
%   (Waarbij C een constante is en 'p' de convergentie-orde).
%
%   Als we aan beide kanten de logaritme nemen, krijgen we:
%       log(res_{k+1}) ≈ log(C) + p * log(res_k)
%
% Dit heeft de vorm van een rechte lijn (y = b + a*x). In een log-log 
% plot is de RICHTINGSCOËFFICIËNT (de helling) van de lijn dus exact 
% gelijk aan de convergentie-orde 'p'.
%
% DEFINITIES VAN CONVERGENTIE:
% - Lineair (p=1):      res_{k+1} ≈ C * res_k^1
%                       De fout neemt per stap met een vaste factor af.
%                       (De lijn in de loglog-plot heeft helling 1).
% - Kwadratisch (p=2):  res_{k+1} ≈ C * res_k^2
%                       Het aantal correcte decimalen verdubbelt ongeveer 
%                       bij elke iteratie-stap zodra we dicht bij de 
%                       oplossing zijn. (Helling in plot is 2).
% - Kubisch (p=3):      res_{k+1} ≈ C * res_k^3
%                       Het aantal correcte decimalen verdrievoudigt.
%                       (Helling in plot is 3).
%
% HOE LEZEN WE DIT AF IN DEZE PLOT?
% 
%   We plotten als referentie drie stippellijnen met bekende hellingen: 
%       t       (helling 1), 
%       t.^2    (helling 2),
%       t.^3    (helling 3).
% 
% Door aan het einde van het convergentieproces (links in de grafiek, 
% waar de fout klein is) te kijken naar de helling van A1, A2 en A3, 
% kunnen we ze simpelweg visueel vergelijken met de stippellijnen:
%
%   - Loopt de curve parallel aan de t stippellijn?    -> Lineair
%   - Loopt de curve parallel aan de t.^2 stippellijn? -> Kwadratisch
%   - Loopt de curve parallel aan de t.^3 stippellijn? -> Kubisch
%%%

% Residus uitzetten TOV vorige residus
loglog(res1(1:n1-1), res1(2:n1), '+-',...
        res2(1:n2-1), res2(2:n2), 'x-',...
        res3(1:n3-1), res3(2:n3), '.-');

hold on

% Twee punten op de x-as
t = [1e-5, 1];

% Referentielijnen voor 3 convergentie-ordes
%   -> log(t^p) = p * log(t)
%   -> rico stelt de convergentie-orde voor
plot(t, t.^3,'--', t, t.^2,'--', t, t, '--')

hold off
legend('A1','A2','A3','Kubisch', 'Kwadratisch', 'Lineair', 'Location','northwest');
xlabel('res_k');
ylabel('res_{k+1}');

%%%
% CONCLUSIE
%
%   A_1 convergeert KUBISCH
%       -> SYMMETRISCH!
%
%   A_2 convergeert KWADRATISCH
%       -> niet symmetrisch
%       -> niet defectief
%
%   A_3 convergeert LINEAIR
%       -> niet symmetrisch
%       -> DEFECTIEF

%%
% V contains the eigenvectors as columns
% D is a diagonal matrix containing the eigenvalues
[V, D] = eig(A3)
% A3 = V D V^-1

% Ziet er op het eerste zicht oke uit, 
% MAAR twee eigenvectoren zijn identiek
abs(V(:,3)) - abs(V(:,4))

% DUS V^{-1} bestaat NIET!!

%%%
% WAAROM CONVERGEERT EEN DEFECTIEVE MATRIX (A3) LINEAIR?
%
%   A3 heeft een Jordan-blok (gecreëerd door L3(3,4) = 1).
%       -> De eigenwaarde 3 komt twee keer voor (algebraïsche 
%          multipliciteit = 2).
%       -> Er is echter maar één onafhankelijke eigenvector voor 
%          deze eigenwaarde (meetkundige multipliciteit = 1).
%       -> Hierdoor is de matrix DEFECTIEF (niet diagonaliseerbaar).
%
%   Gevolg voor het QR-algoritme:
%       -> Het QR-algoritme leunt normaal gesproken op het feit dat 
%          eigenvectoren onafhankelijk zijn. Door slimme 'shifts' toe
%          te passen, scheidt het de eigenvectoren razendsnel (wat zorgt 
%          voor kwadratische convergentie).
%       -> Omdat A3 onafhankelijke eigenvectoren mist, faalt deze 
%          shift-strategie voor de dubbele eigenwaarde.
%       -> Het algoritme vervalt daardoor in wezen naar een standaard 
%          machtsmethode (power iteration). De fout neemt per stap 
%          slechts met een constante factor af.
%       -> Dit resulteert in LINEAIRE convergentie.
%%%