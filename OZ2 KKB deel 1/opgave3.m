%%% Exercise 3

%% Generate data
close all;
x = linspace(-1, 1, 20);
f = exp(x);

% uniforme gewichten
w1 = ones(size(x));

% Chebyshev gewichten
w2 = 1./sqrt(1-x.^2);

% Normaal gezien oneindig, maar hier gewoon 'groot'...
w2(1) = 11; w2(end) = 11;

%% Plot residu for uniform weight function

% Degree of approximation -> VT van graad 5
M = 5;
figure;

%
% kkb1 (overdetermined without normal equations)
% with uniform weight function
%
c = kkb1(x,f,w1,M);

% Evaluate polynomial in all points
y = polyval(c(end:-1:1), x); 

% Compute residu
r = y-f; 

subplot(2,2,1);
plotres(x,r,w1);
title('kkb1')
ylabel('w=1')

%
% kkb2 
%   (efficient version of normal equations where G=(A^T W A) 
%    is computed from the first row and first column of the 
%    Hankel matrix)
% with uniform weight function
%
c = kkb2(x,f,w1,M);

% Evaluate polynomial in all points
y = polyval(c(end:-1:1), x); 

% Compute residu
r = y-f; 

subplot(2,2,2);
plotres(x,r,w1);
title('kkb2')

% 
% kkb1 with Chebychev weight function
%
c = kkb1(x,f,w2,M);
y = polyval(c(end:-1:1), x);
r = y-f;
subplot(2,2,3);
plotres(x,r,w2);
ylabel('w=1/sqrt(1-x^2)')

% 
% kkb2 with Chebychev weight function
%
c = kkb2(x,f,w2,M);
y = polyval(c(end:-1:1), x);
r = y-f;
subplot(2,2,4);
plotres(x,r,w2);

%% Plot max residu for increasing degree

% Graden die gebruikt worden bij benadering
maxGraad = 19;

% Hou residuen bij voor KKB1 en KKB2
residual = zeros(1, M);
residual2 = zeros(1, M);

% Loop over benaderingsgraad van 1 tem maxGraad
for M = 1:maxGraad

    % Overgedetermineerd met uniform gewicht
    c = kkb1(x,f,w1,M);
    y = polyval(c(end:-1:1), x);
    residual(M) = max(abs(y-f));

    % (Hankel) normaal vergelijkingen met uniform gewicht
    c = kkb2(x,f,w1,M);
    y = polyval(c(end:-1:1), x);
    residual2(M) = max(abs(y-f));

end

% Plot
figure;
semilogy(1:maxGraad, residual, 'b.-', 1:maxGraad, residual2, 'r.-');
legend('kkb1', 'kkb2')
grid on;
