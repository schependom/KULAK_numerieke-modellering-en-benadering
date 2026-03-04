%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inleiding Matlab %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Een commentaar wordt geschreven met een % (i.p.v. # in Python).
% Je kan je code opsplitsen in secties en die apart runnen 
% (met 'Run Section'). Het maken van een sectie doe je met 2 %:
%% Basiscommando's
2+2             % Het resultaat wordt altijd geprint,
2+3;            % behalve als je een ; op het einde plaatst.
ans             % Het laatst bekomen resultaat kan je met 'ans' verkrijgen. 
a = 4;          % Assignment is zoals in Python
disp('I <3 Matlab'); % Printen doe je met disp
disp(a);
clear a;        % In je Workspace zie je de opgeslagen variabelen.
                % Die kan je verwijderen door clear te doen.

help clear      % Als je niet weet wat iets doet, typ help ...  

%% Vectoren en matrices
vec1 = [1,5,42,15,6,8]  % Vierkante haakjes!
vec1 = [1 5 42 15 6 8]  % met of zonder comma
mat1 = [1,2,3;4,5,6]    % ; geeft een nieuwe rij aan, wat toelaat matrices te maken.
vec2 = 1:10             % range(1, 10) in Python
vec2 = 1:3:10           % range(1, 10, 3), geeft dus [1, 4, 7, 10]
transpose = mat1.'      % .' is transpose, ' is hermitisch transpose

% Elementen selecteren
element23 = mat1(2,3)   % Ronde haakjes voor indexering!
rij2 = mat1(2,:)        % Merk op dat indexering van 1 begint, niet 0.
kolom3 = mat1(:,3)
selectie = vec1([1:3,5:6])
vec1(1:2) = [3000,200]

% Bewerkingen met matrices en vectoren

A = randi(10,4,3) % 4x3 matrix met random getallen van 1 t.e.m. 10
B = randi(10,4,3)
t = rand(4,1)     % 4x1 matrix met random getallen van 0 tot 1
A+B
A-1
A/2     % elementsgewijs
% A*B     % matrixvermenigvuldiging
A*B.'

x = A\t % Los Ax=t op
format long % Toon meer cijfers
x
format short % Default

% Elementsgewijs
A*A'
A.*A
A.^2
A./A 

% Ingebouwde functies

svd(A) % Returnt singuliere waarden
help svd % Of Matlab website, Google, Mathworks forum, stackoverflow, vriendjes, Raphaël
[U,S,V] = svd(A)
[U,S,V] = svd(A,0) % Gereduceerde SVD (0 == "econ")
[~,S,V] = svd(A,0) % Tilde voor ongewenste uitvoer

vec1
mean(vec1)
mat1
mean(mat1)          % Kolom gemiddelden
mean(mat1, 2)       % Rij gemiddelden
mean(mat1,[1,2])    % Gemiddelde over rijen en kolommen samen
mean(mean(mat1))    % Globaal gemiddelde
mean(mat1(:))       % Globaal gemiddelde (flattening)

C = ones(3)
C = zeros(3,2)
C = eye(3,5)
size(C)

%% Plotten

x = 1:0.05:10;
y = sin(x);
plot(x,y);

figure(); % Nieuwe figuur
plot(x,x.^2);
hold on;  % Meerdere lijnen plotten op dezelfde figuur 
plot(x,x.^3,'m--');
hold off;
title('Mijn mooiste plotje');
xlabel('Tijd');
ylabel('Temperatuur');
legend('Kwadraat', 'Derdemacht');

figure();
semilogy(x,x.^2);
loglog(x,x.^2);

%% Functies

help doeietsnuttigs
[som, prod] = doeietsnuttigs(2,3,4)

% For-lus
s = 0;
for i = 1:10
    s = s+i;
end
s
s = sum(1:10)

% While-lus
t = 0;
while t ~= 10 % Niet /=, niet |=, niet !=
    t = t+1;
end
t

s = 0;
% If-voorwaarde
if s >= 10 % Groter of gelijk aan
    disp('Amai, zo groot!');
elseif s == 0 % Gelijk aan
    disp('NUL!')
elseif s == 1 || s == -1 % || is logische of, && is logische en
    disp('Absolute waarde is 1')
else
    disp('Niets speciaals');
end

%% Debuggen
% Om te debuggen, klik je op de regelnummer waar je wilt dat de executie
% stopt. Doe bijvoorbeeld functie doeietsmindernuttigs open en klik op
% regel nummer 3. Run dan deze sectie.
[som, prod] = doeietsmindernuttigs(2,3,4)

% Timings: time.m


%% Workspace opslaan
save('MijneerstescriptjevoorNMB')
clear all
load('MijneerstescriptjevoorNMB')
% Willekeurige data inladen?
A = importdata('voorbeeld.jpg'); % Kies zelf een afbeelding.
image(A)

