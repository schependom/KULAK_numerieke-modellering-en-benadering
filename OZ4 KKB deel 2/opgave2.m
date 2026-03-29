%%% Exercise 2

%% Plot Inbev data
close all;
figure();
inbev = load('abinbev');
data = inbev.data;
plot(data(:,1), data(:,2));

%% Approximate with polynomial

% Set up x-values, function values, weights and degree
x = data(:,1);          % x-values
f = data(:,2);          % Function values
w = ones(size(x));      % Weights
M = 4;

% Compute approximation
c = kkb1(x,f,w,M); % Compute coefficients 
finter = polyval(c(end:-1:1), x); % Evaluate approximation

%% Evaluate in next days by extrapolation

next = 365;                             % Set how many days to extrapolate
xextra = x(end)+(1:next)';              % Set x-values for extrapolation
fextra = polyval(c(end:-1:1), xextra);  % Evaluate approximation

%% Plot approximation and extrapolation together

hold on;
plot(x, finter, 'g')
plot(xextra, fextra, 'r')