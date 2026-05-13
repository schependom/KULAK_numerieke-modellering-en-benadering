% Gelijktijdige iteratie of ‘simultaneous iteration’ past de 
% methode van de machten toe op meerdere kolommen tegelijk. 
% 
% Op die manier kunnen de d GROOTSTE eigenwaarden en 
% bijhorende eigenvectoren bepaald worden
%
% Schrijf een variant van het algoritme uit om de
% d KLEINSTE eigenwaarden en eigenvectoren van 
% een niet-singuliere matrix A te bepalen.

A = A1;

%% Gewone gelijktijdige iteratie (d GROOTSTE)

Q = orth(rand(4,2));

for k = 1:100
    Z = A * Q;          % vermenigvuldig
    [Q, ~] = qr(Z,0);   % orthogonalisatie
end

lambda = diag(Q'*A*Q)'

%% Aangepaste gelijktijdige iteratie (d KLEINSTE)

Q = orth(rand(4,2));

for k = 1:100
    Z = A \ Q;          % A Z = Q  =>  Z
    [Q, ~] = qr(Z,0);   % orthogonalisatie
end

lambda = diag(Q'*A*Q)'
