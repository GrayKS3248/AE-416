function Vortex_Element_Solver(num_panels)

    panels = get_panels(15.957,0.2025,num_panels);
    control_points = get_control_points(panels);
    circulation_points = get_circulation_points(panels);
    
    visualize(panels, control_points, circulation_points)
    
    
    
end

function panels = get_panels(k,m,num_panels)
    
    x = linspace(0.0,1.0,num_panels+1);
    x1 = x(x <= m);
    x2 = x(x > m);
    y1 = (k/6) * (x1.^3 - 3*m*x1.^2 + m^2*(3-m)*x1);
    y2 = ((k*m^3)/6) * (1 - x2);
    panels = [x;[y1 y2]];
    
end

function control_points = get_control_points(panels)

    delta_x = (circshift(panels(1,:),-1) - panels(1,:));
    delta_y = (circshift(panels(2,:),-1) - panels(2,:));
    x = panels(1,:) + 0.75 * delta_x;
    y = panels(2,:) + 0.75 * delta_y;
    control_points = [x(1:end-1);y(1:end-1)];

end

function circulation_points = get_circulation_points(panels)

    delta_x = (circshift(panels(1,:),-1) - panels(1,:));
    delta_y = (circshift(panels(2,:),-1) - panels(2,:));
    x = panels(1,:) + 0.25 * delta_x;
    y = panels(2,:) + 0.25 * delta_y;
    circulation_points = [x(1:end-1);y(1:end-1)];
    
end

function visualize(panels, control_points, circulation_points)

    plot(panels(1,:),panels(2,:),'color','k');
    hold on
    scatter(panels(1,:),panels(2,:),'filled','o','k')
    scatter(control_points(1,:), control_points(2,:),'+','b');
    scatter(circulation_points(1,:), circulation_points(2,:),'x','r');
    hold off
    xlim([0.0,1.0])
    ylim([-0.04,0.04])
    title("Vortex Element Solver Airfoil")
    xlabel('x/c')
    ylabel('y/c')
    
end