Mcr =0.5:0.00001:0.90;
g1 = -0.6 * (1-Mcr.^2).^-0.5;
g2 = (2*(1.4*Mcr.^2).^-1) .* (((1+(1.4-1)/2)*(1+((1.4-1)/2)*Mcr.^2).^-1).^(1.4/(1-1.4))-1);

plot(Mcr,g1)
hold on
plot(Mcr,g2)
hold off
xlabel("M_{cr}")
ylabel("C_{p,cr}")
title("Critical Mach Number Solution")
legend('Prandlt-Glauert','Isentropic','Location','sw')
grid on
exportgraphics(gcf,'Critical Mach 1.png','Resolution',500)
    
res = abs(g2 - g1);
g1_soln = g1(find(res == min(abs(g2 - g1))));
g2_soln = g2(find(res == min(abs(g2 - g1))));

error = abs(g2_soln - g1_soln)
Mcr = Mcr(find(res == min(abs(g2 - g1))))
cp_crit = g1_soln