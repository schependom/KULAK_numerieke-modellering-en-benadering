%% Opgave 2
%
% Lees een foto in en geef die weer.

f = imread('lena_gray_512.tif');
% f = imread('fourier.tif'); % Set k2 = 200 for this image

N1 = size(f,1);
N2 = size(f,2);

subplot(2,3,1)
imshow(f)
title('Original Image');

% Maak een 2D frequentievoorstelling van de foto 
% en comprimeer de foto op vijf verschillende manieren:
%
% 1. Verwijder de DC componenten (zet op nul)
% 2. Verwijder de k LAAGSTE AC frequenties in beide richtingen (k KLEIN)
% 3. Verwijder de k HOOGSTE AC frequenties in beide richtingen (k GROOT)
% 4. Verwijder de k HOOGSTE frequenties in 1 richting
% 5. Verwijder de k hoogste frequenties in de andere richting.

% De afbeelding kan gecomprimeerd worden door 
% hoge frequenties te verwijderen, 
% of m.a.w., door deze op nul te zetten. 

% Net als bij de 1D FFT zitten de HOGE frequenties in het MIDDEN, 
% rond de Nyquist frequentie. 
% 
% In het geval van een 2D FFT van een (m×n) matrix hebben we dus 
% TWEE BANDEN MET HOGE FREQUENTIES: 
%
%       -> de rijen     rond m/2 
%       -> de kolommen  rond n/2.

%% Compute 2D fft of image

F = fft2(f);

%% 1. Verwijder DC 

Fdc = F;
Fdc(1,1) = 0;

% Plot
subplot(2,3,2)
imshow(ifft2(Fdc)); % Plot inverse fft
title('1. Without DC Component');

%% 2. Verwijder k laagste frequenties (k klein)

k = 3; % Number of frequencies that are filtered out

Fl = F;
Fl([1:k+1, end-k+1:end],:) = 0;
Fl(:, [1:k+1, end-k+1:end]) = 0;

% Plot
subplot(2,3,3)
imshow(uint8(ifft2(Fl))); % Mainly sharp edges remain visible
title('2. Verwijder k laagste frequenties (k klein)')

%% 3. Verwijder k hoogste frequenties (k groot)

k1 = 240; % Number of frequencies that are filtered out (vertical)
k2 = 200; % Number of frequencies that are filtered out (horizontal)

Fh = F;
Fh(floor(N1/2)+1+(-k1:k1),:) = 0;
Fh(:,floor(N2/2)+1+(-k2:k2)) = 0;

% Plot
subplot(2,3,4)
imshow(uint8(real(ifft2(Fh)))); % Sharp edges are smoothed
title('3. Verwijder k hoogste frequenties (k groot) in beide richtingen')

%% 4. Verwijder hoge frequenties in 1 richting

k1 = 240; % Number of frequencies that are filtered out (vertical)
Fx = F;
Fx(floor(N1/2)+1+(-k1:k1),:) = 0; % Compression in vertical dimension

% Plot
subplot(2,3,5)
imshow(uint8(real(ifft2(Fx))));
title('4. Verwijder hoge frequenties in de x-richting')

%% 5. Verwijder hoge frequenties in de andere richting

k2 = 240; % Number of frequencies that are filtered out (horizontal)
Fy = F;
Fy(:,floor(N2/2)+1+(-k2:k2)) = 0; % Compression in horizontal dimension

% Plot
subplot(2,3,6)
imshow(uint8(real(ifft2(Fy))));
title('5. Verwijder hoge frequenties in de y-richting')

%% Look at frequencies that are kept after compression

figure;
subplot(2,3,1)
imshow(f);
title('Original Image');

subplot(2,3,2)
imagesc(log10(abs(Fdc)))
title('Spectrum: No DC Component');

subplot(2,3,3)
imagesc(log10(abs(Fl)))
title('Spectrum: No Low Frequencies');

subplot(2,3,4)
imagesc(log10(abs(Fh)))
title('Spectrum: No High Frequencies');

subplot(2,3,5)
imagesc(log10(abs(Fx)))
title('Spectrum: No Vertical High Freq.');

subplot(2,3,6)
imagesc(log10(abs(Fy)))
title('Spectrum: No Horizontal High Freq.');