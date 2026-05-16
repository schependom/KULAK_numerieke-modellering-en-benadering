function c = train_nn(xi,y)  
    rng(5000);
    c_init = 0.5*randn(23,1);

    options = optimoptions('lsqnonlin','PlotFcn',@optimplotresnorm);
    [c,finalerr] = lsqnonlin(@(c) kostfunctie(c,xi,y), c_init, [], [], options);
    saveas(gcf,strcat('train_nn_convergentie.png'))
end