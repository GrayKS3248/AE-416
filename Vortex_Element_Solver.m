function [gamma_distribution, cl, cm4c] = Vortex_Element_Solver(num_panels, alpha)

    panels = get_panels(15.957,0.2025,num_panels);
    control_points = get_control_points(panels);
    circulation_points = get_circulation_points(panels);
    
    panel_lengths = vecnorm([circshift(panels(1,:),-1); circshift(panels(2,:),-1)] - panels,1);
    panel_lengths = panel_lengths(1:end-1);
    
    [A, panel_normal] = get_A(panels, control_points, circulation_points);
    B = get_B(panel_normal, alpha);
    
    gamma = zeros(length(circulation_points(1,:)), length(alpha(1,:)));
    for curr_alpha = 1:length(alpha(1,:))
        gamma(:,curr_alpha) = A\B(:,curr_alpha);
    end
    
    cl = 2.0 * sum(gamma);
    cmle = -2.0 * cos(deg2rad(alpha)) .* sum(gamma .* vecnorm(circulation_points,1)');
    cm4c = 0.25 * cl + cmle;
    
    gamma_distribution = gamma ./ panel_lengths';
    
    visualize(panels, control_points, circulation_points, cl, cm4c, alpha, gamma_distribution)
    
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

function [A, panel_normal] = get_A(panels,control_points,circulation_points)
    
    A = zeros(length(control_points(1,:)), length(circulation_points(1,:)));
    panel_normal = zeros(2,length(control_points(1,:)));
    
    for curr_point = 1:length(control_points(1,:))
                    
        panel_normal_dirn = panels(:,curr_point+1) - panels(:,curr_point);
        panel_normal_dirn = [panel_normal_dirn(2); -panel_normal_dirn(1)];
        panel_normal_dirn = panel_normal_dirn / norm(panel_normal_dirn);
        panel_normal(:,curr_point) = panel_normal_dirn;
        
        curr_control_point = control_points(:,curr_point);
            
        for curr_vortex = 1:length(circulation_points(1,:))
            
            r_pq = curr_control_point - circulation_points(:,curr_vortex);
            r_pq_norm = norm(r_pq);
           
            induced_v_dirn = [r_pq(2); -r_pq(1)];
            induced_v_dirn = induced_v_dirn / norm(induced_v_dirn);
            
            cos_delta_pq = dot(induced_v_dirn, panel_normal_dirn);
            
            A(curr_point, curr_vortex) = cos_delta_pq / (2 * pi * r_pq_norm);
            
        end
    end

end

function B = get_B(panel_normal, alpha)

    B = zeros(length(panel_normal(1,:)), length(alpha(1,:)));
    alpha_rad = deg2rad(alpha);
    
    for curr_alpha = 1:length(alpha(1,:))
        
        cos_curr_alpha = cos(alpha_rad(curr_alpha));
        sin_curr_alpha = sin(alpha_rad(curr_alpha));
        V_inf_dirn = [cos_curr_alpha -sin_curr_alpha; sin_curr_alpha cos_curr_alpha] * [1;0];
        
        for curr_point = 1:length(panel_normal(1,:))
        
            V_inf_norm = dot(-1.0 * panel_normal(:,curr_point), V_inf_dirn);
            B(curr_point, curr_alpha) = V_inf_norm;
            
        end
        
    end

end

function visualize(panels, control_points, circulation_points, cl, cm4c, alpha, gamma_distribution)

    n = length(panels(1,:)) - 1;

    figure(1)
    plot(panels(1,:),panels(2,:),'color','k');
    hold on
    scatter(panels(1,:),panels(2,:),'filled','o','k')
    scatter(control_points(1,:), control_points(2,:),'+','b');
    scatter(circulation_points(1,:), circulation_points(2,:),'x','r');
    hold off
    xlim([0.0,1.0])
    ylim([-0.04,0.04])
    title_str = strcat("Vortex Element Solver Airfoil n=", num2str(n));
    title(title_str)
    xlabel('x/c')
    ylabel('y/c')
    legend('panel', 'panel vertex', 'control point', 'vortex')
    save_str = strcat("Panels n=", num2str(n), ".png");
    exportgraphics(gcf,save_str,'Resolution',500)
    
    figure(2)
    colororder({'r','b'})
    yyaxis left
    grid on
    plot(alpha,cl,'color','r')
    ylabel('Cl')
    ylim([-1.5,1.5])
    yyaxis right
    plot(alpha,cm4c,'color','b')
    hold on
    line([0 0], [-0.025,0.025],'color','k')
    line([alpha(1),alpha(end)], [0 0],'color','k')
    hold off
    ylabel('Cm4c')
    ylim([-0.025,0.025])
    ax = gca;
    ax.XTick = alpha;
    ax.GridColor = [0.1,0.1,0.1];
    title_str = strcat("Aerodynamic Characteristics n=", num2str(n));
    title(title_str)
    xlabel('alpha [deg]')
    save_str = strcat("Aerodynamic Characteristics n=", num2str(n), ".png");
    exportgraphics(gcf,save_str,'Resolution',500)
    
    close()
    stairs(panels(1,1:end-1), gamma_distribution(:,1)')
    hold on
    for curr_alpha = 6:5:length(alpha(1,:))
        stairs(panels(1,1:end-1), gamma_distribution(:,curr_alpha)')
    end
    hold off
    xlim([0,1])
    ylim([-3.0,3.0])
    xlabel('x/c')
    ylabel('Circulation')
    title_str = strcat("Circulation Distribution n=", num2str(n));
    title(title_str)
    grid on
    legend('alpha = -10', 'alpha = -5', 'alpha = 0', 'alpha = 5', 'alpha = 10')
    save_str = strcat("Circulation Distribution n=", num2str(n), ".png");
    exportgraphics(gcf,save_str,'Resolution',500)
    
end