% Grid plot
gridN = 500;
Dx1 = 1/gridN;
Dx2 = 1/gridN;
x1vals = [0:Dx1:1];
x2vals = [0:Dx2:1];
for k1 = 1:gridN+1
    for k2 = 1:gridN+1
        x = [x1vals(k1);x2vals(k2)];
        Fx = neuraal_net(x,c);
        Aval(k2,k1) = Fx(1);
        Bval(k2,k1) = Fx(2);
     end
end

[X,Y] = meshgrid(x1vals,x2vals);

p2 = figure;
clf
a2 = subplot(1,1,1);
Mval = Aval>Bval;
contourf(X,Y,Mval,[0.5 0.5])
hold on
colormap([1 1 1; 0.8 0.8 0.8])
plot(xi(1,1:5),xi(2,1:5),'ro','MarkerSize',12,'LineWidth',4)
plot(xi(1,5+1:end),xi(2,5+1:end),'bx','MarkerSize',12,'LineWidth',4)
a2.XTick = [0 1];
a2.YTick = [0 1];
a2.FontWeight = 'Bold';
a2.FontSize = 16;
xlim([0,1])
ylim([0,1])
saveas(p2,strcat('neuraal_net_resultaat.png'))
