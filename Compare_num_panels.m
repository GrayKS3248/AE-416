colororder({'r','b'})
yyaxis left
grid on
plot([0 0],[0 0],'color','k','LineStyle','-')
hold on
plot([0 0],[0 0],'color','k','LineStyle','--')
plot([0 0],[0 0],'color','k','LineStyle','-.')
plot([0 0],[0 0],'color','k','LineStyle',':')
plot(alpha,cl_20,'color','r','LineStyle','-','MarkerSize',0.0001)
plot(alpha,cl_40,'color','r','LineStyle','--','MarkerSize',0.0001)
plot(alpha,cl_60,'color','r','LineStyle','-.','MarkerSize',0.0001)
plot(alpha,cl_80,'color','r','LineStyle',':','MarkerSize',0.0001)
hold off
ylabel('Cl')
ylim([-1.5,1.5])
yyaxis right
plot(alpha,cm4c_20,'color','b','LineStyle','-')
hold on
plot(alpha,cm4c_40,'color','b','LineStyle','--')
plot(alpha,cm4c_60,'color','b','LineStyle','-.')
plot(alpha,cm4c_80,'color','b','LineStyle',':')
line([0 0], [-0.025,0.025],'color','k')
line([alpha(1),alpha(end)], [0 0],'color','k')
hold off
ylabel('Cm4c')
ylim([-0.025,0.025])
ax = gca;
ax.XTick = alpha;
ax.GridColor = [0.1,0.1,0.1];
title("Aerodynamic Characteristics")
legend('n=20', 'n=40', 'n=60', 'n=80','Location','nw')
xlabel('alpha [deg]')
exportgraphics(gcf,'Aerodynamic Characteristics.png','Resolution',500)

close()
