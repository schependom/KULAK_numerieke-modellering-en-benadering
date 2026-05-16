function [H, Q, rw] = arnoldiritz(A, b, N)
% function [H, Q, rw] = arnoldiritz(A, b, N)
%
%%%%% BEGIN 
%%%%% ZELFDE ALS ArnoldiEW.m
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
%%%%% EINDE
%%%%% ZELFDE ALS ArnoldiEW.m
%
%   rw      Matrix van dimensie (N x N) die de evolutie van alle 
%           Ritz-waarden bijhoudt. 
%           - Kolom n bevat de n berekende Ritz-waarden van stap n.
%           - De overige rijen in die kolom (n+1 tot N) zijn NaN.

% Preallocation
% (We initialiseren rw met NaN zodat lege plekken in de matrix geen 
%  nul-waarden worden, wat verwarrend kan zijn bij het plotten).
rw = nan(N,N); 

m = size(A,1);     
Q = zeros(m,N+1);
H = zeros(N+1,N);

% Initialisatie
Q(:,1) = b(:)/norm(b);

for n = 1:N
    
    v = A*Q(:,n);
    for j = 1:n
        H(j,n) = Q(:,j)'*v;
        v = v - H(j,n)*Q(:,j);
    end
    H(n+1,n) = norm(v);
    
    if H(n+1,n) <= 0
        break; 
    end
    
    Q(:,n+1) = v/H(n+1,n);
    
    % Bereken ALLE Ritz-waarden van de huidige H_n
    ewH = eig(full(H(1:n,1:n)));
    
    % Sla deze n Ritz-waarden op in de n-de kolom van de rw matrix.
    % (We berekenen hier GEEN fout t.o.v. de echte eigenwaarden van A)
    rw(1:n,n) = ewH;
    
end