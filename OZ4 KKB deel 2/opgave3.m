%%% Exercise 3

%% Initialize stuff

Mmax = 20; % Maximal approximation degree
res = zeros(Mmax, 4); % Matrix holding all norms of residuals

%% Plot arcsine

x = linspace(0,1,100)';
f = asin(x);
close all;
figure();
plot(x, f);

%% First approximation: 
% equidistant points, w = 1

ex = 1;
x = linspace(0,1,100)';
f = asin(x);
w = ones(size(x));
approx = zeros(length(x),Mmax);

% Compute approximation for different degrees
for M = 1:Mmax
    c = kkb1(x,f,w,M);
    approx(:,M) = polyval(c(end:-1:1),x);
    res(M,ex) = norm(f-approx(:,M));
end

% Plot degree five approximation
hold on;
plot(x,approx(:,5));

%% Second approximation: 
% equidistant points, Chebychev weight function

ex = 2;
x = linspace(0,0.9999,100)';
f = asin(x);
w = 1./sqrt(1-x.^2);
approx = zeros(length(x),Mmax);

% Compute approximation for different degrees
for M = 1:Mmax
    c = kkb1(x,f,w,M);
    approx(:,M) = polyval(c(end:-1:1),x);
    res(M,ex) = norm(f-approx(:,M));
end

% Plot degree five approximation
plot(x,approx(:,5));

%% Third approxmation: 
% cosine points, w = 1

ex = 3;
x = cos(linspace(-pi/2,0,100))';
f = asin(x);
w = ones(size(x));
approx = zeros(length(x),Mmax);

% Compute approximation for different degrees
for M = 1:Mmax
    c = kkb1(x,f,w,M);
    approx(:,M) = polyval(c(end:-1:1),x);
    res(M,ex) = norm(f-approx(:,M));
end

% Plot degree five approximation
plot(x,approx(:,5));
legend('arcsin','Standard', 'Chebychev', 'Other points'); 

%% Plot residuals for different degrees

figure();
semilogy(1:Mmax,res(:,1));
hold on;
semilogy(1:Mmax,res(:,2));
semilogy(1:Mmax,res(:,3));

legend('Standard', 'Chebychev', 'Other points'); 

%% Compute approximation of special form

ex = 4;
x = linspace(0,0.9999,100)';
w = ones(size(x)); % Weights
f = asin(x);

% Approximate transformation of f with polynomial
ftrans = (pi/2 - asin(x))./sqrt(1-x); % Isolate polynomial part of approximation

for M = 1:Mmax

    % Compute polynomial approximation of transformed f
    c = kkb1(x,ftrans,w,M); 

    % Evaluate polynomial part of approximation
    polyapprox = polyval(c(end:-1:1),x); 

     % Transform again to approximation of f
    approx(:,M) = -sqrt(1-x).*polyapprox+pi/2;

    % Compute residual
    res(M,ex) = norm(f-approx(:,M));

end

%% Plot approximation of special form

% Plot degree five approximation 
figure(1);
hold on;
plot(x, approx(:,5));
legend('arcsin','Standard', 'Chebychev', 'Other points','Special'); 

% Plot residual
figure(2);
hold on;
semilogy(1:Mmax, res(:,4));
legend('Standard', 'Chebychev', 'Other points','Special'); 