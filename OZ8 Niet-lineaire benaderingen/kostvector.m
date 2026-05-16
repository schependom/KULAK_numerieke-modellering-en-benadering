function kostvector = kostfunctie(c,xi,y)
    yi = neuraal_net(xi, c);
    kostvector = 1/2 * vecnorm(y-yi, 2, 1);
end