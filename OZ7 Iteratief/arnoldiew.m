function [H, Q, maxd, mind] = arnoldiew(A, b, N)
% function [H, Q, maxd, mind] = arnoldiew(A, b, N)
%
% INPUTS
%   
%   A       vierkante matrix van dimensie (m x m)
%   b       een willekeurige startvector van dimensie (m x 1)
%   N       het aantal iteratiestappen dat we willen uitvoeren (N <= m)
%
% OUTPUTS
%
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


% De functie full(A) pakt een (mogelijk) sparse matrix en 
% converteert deze naar een 'normale' volle matrix in het geheugen.
%
% eig() werkt niet met spaarse matrices, dus we moeten full() aanroepen
% om hem volledig in het geheugen op te slaan


% Echte eigenwaarden van volle matrix berekenen ter referentie
ewA = eig(full(A));
maxewA = max(ewA);
minewA = min(ewA);

% Preallocation
m = size(A,1);      % aantal rijen/kolommen van A (dimensie m)
Q = zeros(m,N+1);   % matrix voor de N+1 orthonormale basisvectoren
H = zeros(N+1,N);   % uitgebreide boven-Hessenberg matrix ~H_N

% Initialisatie (normaliseer de startvector)
Q(:,1) = b(:)/norm(b); % q_1

for n = 1:N
    
    v = A*Q(:,n);
    
    % Gemodificeerd Gram-Schmidt proces voor orthogonalisatie
    for j = 1:n
        H(j,n) = Q(:,j)'*v;
        v = v - H(j,n)*Q(:,j);
    end
    
    H(n+1,n) = norm(v);
    
    % Check voor 'lucky breakdown'
    %   -> gevonden Krylov-deelruimte is zelfde als vorige
    %   -> Ritz-waarden (ew van H_n) zijn geen benaderingen meer, maar
    %      zijn de EXACTE eigenwaarden van de matrix A.
    if H(n+1,n) <= 0
        break; 
    end
    
    % Volgende Arnoldi vector
    Q(:,n+1) = v/H(n+1,n);
    
    % Bereken Ritz-waarden van H_n (de vierkante n x n submatrix)
    ewH = eig(full(H(1:n,1:n)));
    maxewH = max(ewH);
    minewH = min(ewH);
    
    % Foutregistratie
    maxd(n) = abs(maxewA - maxewH);
    mind(n) = abs(minewA - minewH);
    
end