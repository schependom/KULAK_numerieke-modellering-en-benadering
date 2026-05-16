function Fx = neuraal_net(x,c)
    [W2,W3,W4,b2,b3,b4] = lees_parameters(c);
    sigma2 = activatiefunctie(x, W2, b2);
    sigma3 = activatiefunctie(sigma2, W3, b3);
    Fx = activatiefunctie(sigma3, W4, b4);
end