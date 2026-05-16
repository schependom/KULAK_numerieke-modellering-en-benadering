% Implementeer de Arnoldi iteratie. 
% 
% De iteratie hoeft geen rekening te houden met ‘breakdown’ (deling door 0). 
% 
% Bepaal bij elke iteratie de grootste en kleinste Ritz-eigenwaarden 
% (efficientie is niet belangrijk, gebruik de standaard Matlab-functie). 
% 
% Construeer een matrix met 10 eigenwaarden in het interval (4,5). 
% Maak een grafiek van de absolute waarde van het ver- schil tussen de 
% grootste Ritz-eigenwaarde en de grootste eigenwaarde. Idem voor de 
% kleinste Ritz-eigenwaarde en de kleinste eigenwaarde. 
% 
% Doe nu hetzelfde maar vervang een van de eigenwaarden door 8 
% en vervolgens ook een door 2. Wat stel je vast? 
% Wat zou je hieruit kunnen besluiten?

%%
m = 10;
L = eigint(4,5,m);
L(1) = 8;
L(2) = 2;


% willglv(L):
%   pas een willekeurige gelijkvormigheidstransformatie toe op 
%   een matrix met als diagonaalelementen L(:)

[A,V] = willglv(L);

%%%%%%%%%%%%
%%% ARNOLDI
%%%%%%%%%%%%
%
%
% !! IMPORTANT !!
%
%   -> A is (m x m)
%   -> b is (m x 1)
%
% For the rest of the book!!
%
%
% COMPLETE reduction to Hessenberg:
%
%   A = Q H Q^*     <=>     AQ = QH
%
%
% If m (#rows) is HUGE or infinite, this is not possible.
% 
%   => We consider the FIRST n COLUMNS (!) of AQ=QH
%
%   => Q_n = (m x n) matrix with first n columns of Q
%
%   => ~H_n = (n+1 x n) upper-left section of H
%           = also Hessenberg
%
%
% So, we have:
%
%   => A Q_n = Q_{n+1} ~H_n
%
%   => A q_n = (h_1n q_1 + ... + h_nn q_n) + h_{n+1,n} q_{n+1}
%
%   
% So we see that q_{n+1} satisfies an (n+1)-term recurrence relation
% involving itself and the previous Krylov vectors q_1, ..., q_n


%%%
% STEP 1 in Arnoldi Algorithm (p 252)
%
%   1. b is an arbitrary vector
%
%   2. we normalize b:
%
%           q_1 = b/||b||

b = rand(m,1); % random m-vector met waarden (0,1)


%   H       De uitgebreide boven-Hessenberg matrix ~H_N van dimensie (N+1 x N).
%             - H(1:N, 1:N) is de projectie H_N = Q_N^* A Q_N op de Krylov-deelruimte.
%             - De matrix voldoet aan de relatie: A * Q_N = Q_{N+1} * ~H_N
%             - Als N = m, dan bevat H(1:m, 1:m) de volledige Hessenberg-vorm van A.
%
%   Q       Matrix van dimensie (m x N+1) met orthonormale kolommen.
%           De eerste N kolommen (Q_N) spannen de Krylov-deelruimte K_N(A,b) op.
%           De (N+1)-de kolom is nodig voor de residu-berekening in de volgende stap.
%
%   maxd(:) array (lengte N) met de absolute waarde van het verschil tussen de
%           echte grootste eigenwaarde van A en de grootste Ritz-waarde in stap n.
%   mind(:) array (lengte N) met de absolute waarde van het verschil tussen de
%           echte kleinste eigenwaarde van A en de kleinste Ritz-waarde in stap n.
[H, Q, maxd, mind] = arnoldiew(A, b, m);


%%
t = 1:m;
semilogy(t, maxd, '+-', t, mind, 'x-');
xlabel('Iteratie');
ylabel('Absolute fout');
legend('Max eigenwaarde','Min eigenwaarde')

% Save plot as PDF
% exportgraphics(gcf, 'eigenwaarden-in-45-en-8-en-2.pdf', 'ContentType', 'vector');

%%% 
% BESLUIT
%
% Extreme eigenwaarden worden sneller goed benaderd!!