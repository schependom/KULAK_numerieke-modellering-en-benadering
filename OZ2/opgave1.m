%%% Exercise 1

%% Initialize 20 random points, function values and weights
x = linspace(-1, 1, 20);
r = randn(size(x));
w = rand(size(x));

% Run plotres function
close all;
figure;
box on
plotres(x,r,w)
