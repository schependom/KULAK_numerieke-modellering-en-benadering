%%% Exercise 1a

%% Plot f
x = linspace(0.0001,2,100);
f = expint(1) - expint(x);
close all;
figure();
plot(x,f); % gaat altijd maar dichter naar 0