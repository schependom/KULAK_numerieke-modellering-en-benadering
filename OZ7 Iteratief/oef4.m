%%%
%
%
% Iteratieve methodes zijn in het bijzonder geschikt wanneer de 
% matrix een IJLE structuur heeft.
%
%
% Genereer een random, ijle matrix (1000×1000) met het commando sprand.
% Pas een aantal Arnoldi-iteratiestappen toe op deze matrix (bv. 100). 
% 
% Bereken per iteratiestap de Ritz waarden. Daarbij mag je gebruik maken 
% van het ingebouwde Matlab-commando grafiek waarin je toont hoe de 
% Ritz waarden per iteratiestap convergeren naar de eigenwaarden van A. 
% 
% Toon daarbij enkel de reele delen. 
% Bespreek bondig het convergentiegedrag.
%
%
%%%



% Aantal uitgevoerde Arnoldi iteraties
N = 100;

% SParse RANDom matrix (uniformly distributed)
A = sprand(1000,1000, 0.01);

% (m x 1) random beginvector
b = randn(size(A,2),1);

% Eigenwaarden van een SPARSE matrix 
%
%   -> met eigs()
%   -> want eig() werkt alleen met full()
ewA = eigs(A)
%   4.9587 + 0.0000i
%   0.1730 + 1.8258i
%   0.1730 - 1.8258i
%  -1.5798 + 0.9259i
%  -1.5798 - 0.9259i
%  -0.3624 + 1.7809i



% ~H_N, Q_N, ritzwaarden
[H, Q, rw] = arnoldiritz(A,b,N);

% Pak de reële delen en sorteer ze van groot naar klein
gesorteerde_ew = sort(real(ewA), 'descend');

% Wijs de grootste en op-een-na-grootste toe
eig_1 = gesorteerde_ew(1);
eig_2 = gesorteerde_ew(2);

% Bereken alpha
alpha = abs(eig_2) / abs(eig_1);


%%%
% ALPHA bepaalt de convergentiesnelheid van de machtsmethode.
%
% \alpha = | lambda_2 / lambda_1 |
%
% Beschouw de ontbinding
%
%   A^n b = lambda_1^n [c_1 v_1 + c_2 (lambda_2 / lambda_1)^n v_2
%                               + c_3 (lambda_3 / lambda_1)^n v_3
%                               + ...
%                               + c_m (lambda_m / lambda_1)^n v_m
%
% De term die het langzaamst naar nul gaat, 
% is de term met de op-een-na-grootste eigenwaarde.



% Maximale ritzwaarde per iteratiestap 
rw_max = real(max(rw));

% Plotten
figure
% Fout ten opzichte van de ECHTE grootste eigenwaarde (eig_1)
semilogy(abs(rw_max - eig_1), 'LineWidth', 1.5)
hold all

% Theoretische convergentielijn
n = 1:N;
semilogy(n, alpha.^n, '--', 'LineWidth', 1.5);

xlabel('Iteratiestap');
ylabel('Absolute fout');
legend('Fout grootste Ritz-waarde', 'Theoretische grens (\alpha^n)')

%%%
% BESLUIT & CONVERGENTIEGEDRAG
%
% 1. Fase van Lineaire Convergentie (Iteratie 1 t/m ~35)
%    De fout van de grootste Ritz-waarde neemt af met een vaste factor 
%    per iteratiestap (e_{n+1} ~ \alpha * e_n). Binnen de numerieke 
%    wiskunde heet dit 'lineaire convergentie'. Omdat we de fout 
%    plotten op een logaritmische y-as (semilogy), vertaalt deze 
%    lineaire convergentie zich visueel in een strak dalende rechte 
%    lijn. De werkelijke fout daalt minstens even snel als de 
%    theoretische grens (\alpha^n).
%
% 2. Fase van Stagnatie (Iteratie ~35 t/m 100)
%    Rond iteratie 35 stuit de fout op een harde grens van ~10^-15. 
%    Dit is de MACHINEPRECISIE van Matlab (double-precision). Het 
%    verschil tussen de eigenwaarde en de Ritz-waarde is letterlijk 
%    te klein geworden voor de computer om nog te registreren. Verdere 
%    iteraties hebben hier geen zin meer en leiden enkel tot meer 
%    verlies van orthogonaliteit.
%
%
% DE LINK TUSSEN ARNOLDI EN DE MACHTSMETHODE
%
%   De theoretische grens (\alpha^n) met \alpha = |lambda_2 / lambda_1|
%   is eigenlijk de convergentiesnelheid van de MACHTSMETHODE.
%
%   Verschil in aanpak:
%   -> De machtsmethode berekent de iteratieslag b, Ab, A^2b... en kijkt 
%      uiteindelijk ENKEL naar de allerlaatste vector A^{n-1}b. Alle 
%      informatie uit de eerdere stappen wordt weggegooid.
%
%   -> Arnoldi bouwt precies dezelfde reeks op, maar BEWAART alle 
%      voorgaande vectoren als basis voor de Krylov-deelruimte. Arnoldi 
%      zoekt vervolgens in die HELE opgespannen ruimte naar de optimale 
%      benadering (de Ritz-waarden).
%
%   Gevolg voor convergentie:
%
%   Omdat Arnoldi optimaliseert over de gehele n-dimensionale Krylov 
%   ruimte, in plaats van slechts een 1-dimensionale vector, fungeert 
%   de machtsmethode als de 'worst-case' THEORETISCHE BOVENGRENS.
%
%   Arnoldi zal altijd minstens even snel, maar in de praktijk vaak 
%   nóg sneller convergeren.
%%%