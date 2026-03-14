format short e;

%% Matrix
X = rand(5,2)

% We maken de sub-diagonaal elementen van de eerste kolom
% van X (elementen van x_1 onder de diagonaal) extreem klein.
%
%       => ||x_1|| = \approx x_11
%       => de norm van de eerste kolom is ongeveer x_11

X(2:5,1) = X(2:5,1)*1.e-8;

%% Beide reflecties toepassen op X
% De bedoeling is dat de eerste kolom van X alleen maar
% nullen bevat onder de diagonaal, omdat we spiegelen rond
% het hypervlak dat loodrecht staat op x_1 +/- ||x_1||e_1

% Reflectievectoren v
U1  = X(:,1) + norm(X(:,1))*eye(size(X,1),1)
U2  = X(:,1) - norm(X(:,1))*eye(size(X,1),1) % nul op diagonaal!

H1X = X  - 2 / (U1'*U1) * U1 * U1' * X 
%   -> quasi perfect
H2X = X  - 2 / (U2'*U2) * U2 * U2' * X 
%   -> niet goed!
%   -> omdat X(1,1) en norm(X(:,1)) nagenoeg gelijk zijn. Het aftrekken
%      van deze getallen in U2(1) veroorzaakt catastrofale uitdoving,
%      waardoor dit eerste element verandert in effectief nul (numerieke ruis).
%   -> Hierdoor gedraagt de vector U2 zich intern als [0; subdiagonaal].
%   -> De daaropvolgende Householder-transformatie (H2X) faalt daardoor
%      in zijn doel: in plaats van de subdiagonaal op nul te projecteren,
%      klapt de spiegeling door deze verminkte U2 enkel het teken van de
%      subdiagonaal-elementen om (de kleine foutwaarde blijft dus behouden).

%% We kunnen aantonen dat U2 = UU2 (niet numeriek)
UU2 = -sum(X(2:5,1).^2)/U1(1);
UU2(2:5,1) = X(2:5,1);
H3X = X  - 2 / (UU2'*UU2) * UU2 * UU2'* X 

% Je berekent met UU2 wiskundig exact hetzelfde als U2. 
% Omdat je nu echter kwadraten optelt en deelt door een som van 
% twee positieve getallen (in plaats van ze van elkaar af te trekken), 
% verlies je geen precisie. De foutwaarde in je variabele E zal 
% hierdoor aantonen dat het verschil tussen U2 en UU2 uitsluitend 
% uit deze afrondingsfout bestaat.

E = (U2-UU2)./UU2