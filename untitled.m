figure
x = 0:0.01:4;
V = (1-exp(-5*(x-1.5))).^2;
V(x<=1.3)=nan;
Vc = 2./x+1.5;

plot(x,V,"Color",'k',"LineWidth",2.5)
hold on
plot(x,Vc,"Color",'k',"LineWidth",2.5)
ylim([-0.25,3.75])
xlim([1,4])