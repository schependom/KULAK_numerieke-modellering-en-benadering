%%% 
% INVERSE ITERATIE
%
% x = b / ||b||
%
% for k=1,2, ...
%
%   (A - \mu I) y = x               % y = (A - \mu I)^{-1} x
%               x = y / ||y||


% 1.    x \approx EV            <==         \mu \approx \lambda
%
% 2.    \mu \approx \lambda     ==>         stelsel bijna singulier!
%                                           grote relatieve fouten op y
%
% Toch blijkt de berekende x steeds een goede benadering voor de
% eigenvector in geval (2). Verklaar de schijbare tegenspraak.


% A = A_1
L1 = diag(1:4);
P1 = orth(rand(4));
A1 = P1 * L1 * P1';
A = A1;

eig(A); % 1, 4, 3, 2

% mu
mu = 2 + 1e-5;

% k(A)
cond(A - mu*eye(size(A))) % 2 x 10^5

% Voor stelsel Ax=b:
%
%       ||dx||/||x|| <= k(A) ||db|| / ||b||
%
% Voor relatieve fout van e_mach=10e-16, is de maximale fout op x ongeveer
% 10e-11, wat niet zo groot is.
%
% Each decade of κ costs exactly one decimal digit of accuracy.
%
%       With κ ≈ 2×10⁵, you lose log₁₀(2×10⁵) ≈ 5.3 digits, 
%       leaving ~10–11 correct digits.


% 1 ieratie van het algoritme
Am = A - mu*eye(size(A));
b = rand(4, 1);
x = b/norm(b);
y = Am\x        % Am y = x oplossen ('inverse')
x = y/norm(y)

% We kunnen dit vergelijken met de echte eigenvector 
% behorende bij de eigenwaarde 2.
%
% EVD: A1 = P1 * L1 * P1'
% DUS: P1 bevat de eigenvectoren in de kolommen
%
eigvec = P1(:, 2);
format SHORTG 
[norm(x - eigvec) norm(x + eigvec)] % [2 9e-5]
%
% VERKLARING SCHIJNBARE TEGENSPRAAK
%
% Schrijf de startvector x uit in de eigenbasis {v1, v2, v3, v4} van A:
%
%       x = α1 v1 + α2 v2 + α3 v3 + α4 v4
%
% Dan geldt voor de oplossing y = (A − µI)^{-1} x:
%
%       y =  α1/(λ1−µ) v1  +  α2/(λ2−µ) v2  +  α3/(λ3−µ) v3  +  α4/(λ4−µ) v4
%
% De coëfficiënten 1/(λi−µ) zijn:
%
%       1/(1−µ) ≈ −1       1/(2−µ) ≈ −10^5      1/(3−µ) ≈ 1      1/(4−µ) ≈ 0.5
%
% De component langs v2 wordt dus met een factor ~10^5 uitvergroot,
% terwijl alle andere componenten van orde 1 blijven.
% Na normalisatie domineert de v2-component volledig:
%
%       x = y/||y||  ≈  ±v2
%
% De grote fouten op y liggen dus wél langs v2, maar ook y zelf
% ligt bijna volledig langs v2. De richting van y is robuust,
% ook al is de grootte sterk verstoord.

% Stel we perturberen A-µI:
xp = b/norm(b);

% Am is nog steeds shifted A
Ap = Am + (rand(4)*2 - 1) * 1e-8;
yp = Ap \ xp
xp = yp / norm(yp)

% Wat is de fout?
norm(y-yp) % 9.25
norm(x-xp) % 7e-8
%
% VERKLARING GROTE FOUT OP y, KLEINE FOUT OP x
%
% De perturbatie δA (met ||δA|| ≈ 10^{-8}) veroorzaakt via de slechte
% conditionering (κ ≈ 10^5) een fout op y van orde:
%
%       ||δy|| / ||y||  ≲  κ · ||δA|| / ||A−µI||  ≈  10^5 · 10^{-8}  ≈  10^{-3}
%
% Omdat ||y|| zelf al van orde 10^5 is (door de versterking langs v2),
% is de absolute fout ||δy|| ≈ 10^{-3} · 10^5 = O(1), wat we waarnemen.
%
% Cruciaal: de perturbatie δy heeft, net als y zelf, zijn grootste component
% langs v2 (de richting die (A−µI)^{-1} het sterkst uitvergroot).
% Zowel y als yp wijzen dus nagenoeg in dezelfde richting.
% Na normalisatie valt de gemeenschappelijke v2-component weg uit de fout:
%
%       ||x − xp||  =  ||y/||y|| − yp/||yp||||  ≈  O(10^{-8})
%%%
% BESLUIT
% We zien dat y sterk veranderd is (fout van orde 1),
% terwijl x goed blijft (fout van orde 10−9).
%
% We kunnen besluiten dat de oplossing van het geperturbeerde stelsel
% yp inderdaad sterk verandert, wat te verwachten is aangezien het stelsel
% slecht geconditioneerd is.
%
% De richting van de y en yp blijft echter nagenoeg gelijk,
% wat we te zien krijgen als we de oplossing normaliseren (x).
%
% De fout en de perturbaties worden immers vooral versterkt in de richting
% die we zoeken en verzwakt in de andere richtingen. (Door matrix A− µI.)