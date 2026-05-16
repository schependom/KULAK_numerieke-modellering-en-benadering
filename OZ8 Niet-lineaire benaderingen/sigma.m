function sigma = activatiefunctie(x,W,b)
    sigma = 1./(1+exp(-(W*x+b)));
end