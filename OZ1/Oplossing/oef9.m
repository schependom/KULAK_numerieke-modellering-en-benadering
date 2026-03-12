%% Setup.
close all;
clear all;
clc;

%% Genereer data.
t = linspace(0, 50, 200);
V{1} = [sin(2 * pi * t / 50).', sin(2 * pi * t / 25).'];
V{2} = V{1};
A = V{1} * V{2}.'; % Genereer matrix van factor matrices.

% Bekijk het signaal zonder ruis.
figure(Name = 'Ruisloze data');
surf(A);
title('Ruisloze matrix');

%% Voeg ruis toe.
SNR = 10;          % (signaal-ruisverhouding in decibel)
ruis = randn(size(A));
factor = 10^(-SNR/20) * norm(A, 'fro') / norm(ruis, 'fro');
Ar = A + factor * ruis;

% Bekijk ruizig signaal.
figure(Name = 'Ruizige data');
surf(Ar);
title('Ruizige matrix');

%% Bepaal de benaderingsrang.
% We merken dat er twee singuliere waarden significant zijn,
% dus voeren we een rang-2 benadering uit van de matrix.
[UAr, SAr, VAr] = svd(Ar);
    fprintf('singular waarden (10 grootsten): ');
    fprintf('  %.2f', diag(SAr(1:10,1:10)));
    

%% Voer een lage-rang benadering uit.
% De ruis kan van de data gehaald worden door een
% singulierenwaardenontbinding te berekenen en het resultaat af te kappen.
% Het idee is dat de kleine singuliere waarden die verwacht worden nul te
% zijn op nul worden gezet.

% Vorm lage-rang benadering.
Artrunc = UAr(:,1:2) * SAr(1:2,1:2) * VAr(:,1:2)';

% Vergelijk het resultaat grafisch.
figure(Name = 'Data na denoising');
surf(Artrunc);
title('Matrix met ruis verwijderd');

figure(Name = 'Doorsnede');
title('Ruisloze, ruizig en denoised doorsnede');
plot(A(:,40))
hold
plot(Ar(:,40),'r')
plot(Artrunc(:,40),'g')
legend({'Ruislooos', 'Ruizig', 'Lage-rang benadering'})

%% Frobeniusnorm
function f = frob(A)
    % Berekent de Frobenius-norm van matrix A
    % Dit is hetzelfde als norm(A, 'fro')
    f = norm(A, 'fro');
end


%% Bereken de relatieve nauwkeurigheid.
% Het verschil tussen de ruizige matrix en de lage-rang benadering komt
% overeen met de signaal-ruisverhouding.
esterrAr = 20 * log10(frob(Ar - Artrunc) / frob(Ar));
fprintf('\nRelatief verschil ruizig vs denoised matrix: %.2fdB\n', esterrAr);
% Het relatief verschil vergeleken met de ruizige tensor is ongeveer -10dB,
% wat overeenkomt met de gegeven signaal-ruisverhouding (snr = 10).

esterrA = 20 * log10(frob(A - Artrunc) / frob(A));
fprintf('Relatief verschil ruisloos vs denoised matrix: %.2fdB\n', esterrA);
% Het relatieve verschil tussen het ruisloze matrix en de lage-rang
% benadering van de ruizige matrix is kleiner dan -10dB, wat betekent dat
% de singulierewaardenontbinding een goede benadering geeft van de
% originele matrix.

%%
% Alternatieve berekening, in de wetenschap dat
%
%       v = rang van onze benadering
%       r = rang van A^T A = rang van A
%
%       ||A-A_v||_2 = \sigma_{v+1}
%       ||A-A_v||_F = \sqrt(\sigma_{v+1}^2 + ... \sigma_{r}^2)

% Haal alle singuliere waarden op
sigmas = diag(SAr);

% Bereken de norm van de fout (vanaf index 3 tot het einde)
fout_norm = sqrt(sum(sigmas(3:end).^2));

% Bereken de norm van de volledige matrix Ar
totaal_norm = sqrt(sum(sigmas.^2)); % Of gewoon norm(Ar, 'fro')

% Alternatieve dB berekening
esterrAr_alt = 20 * log10(fout_norm / totaal_norm)

% Verschil
abs(esterrAr_alt-esterrAr)