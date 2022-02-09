function plot_robot(x_ev)

    n_step = 50;
    stop_time = x_ev.time(end);%6;%1.3;
    time_vect = linspace(0,stop_time,n_step);
    
    x_ev_res = resample(x_ev,time_vect);
    
    s_vect = 0:0.0125:1;
    for i_t = 1:length(time_vect)
        xy_c = nan(length(s_vect),2);
        for i_s = 1:length(s_vect)
            xy_c(i_s,:) = forward_kinematics(x_ev_res.Data(i_t,1),x_ev_res.Data(i_t,2),0,s_vect(i_s));
        end
        hold on, plot(-xy_c(:,2),xy_c(:,1), 'linewidth', 0.4, 'color', [230,230,230]/255)
    end
    
    
    for i_s = 1:length(s_vect)
        xy_c_i(i_s,:) = forward_kinematics(x_ev.Data(1,1),x_ev.Data(1,2),0,s_vect(i_s));
        xy_c_e(i_s,:) = forward_kinematics(x_ev.Data(end,1),x_ev.Data(end,2),0,s_vect(i_s));
        xy_c_e_p(i_s,:) = forward_kinematics(x_ev.Data(end,1),x_ev.Data(end,2),+0.5,s_vect(i_s));
        xy_c_e_m(i_s,:) = forward_kinematics(x_ev.Data(end,1),x_ev.Data(end,2),-0.5,s_vect(i_s));
    end
    hold on, plot(-xy_c_e(:,2),xy_c_e(:,1), 'linewidth', 2, 'color', [0, 0.4470, 0.7410])
%     hold on, plot(-xy_c_e_p(:,2),xy_c_e_p(:,1), '--', 'linewidth', 1, 'color', [0, 0.4470, 0.7410])
%     hold on, plot(-xy_c_e_m(:,2),xy_c_e_m(:,1), '--', 'linewidth', 1, 'color', [0, 0.4470, 0.7410])
    hold on, plot(-xy_c_i(:,2),xy_c_i(:,1), 'linewidth', 2, 'color', [0.47,0.67,0.19])
    
    %x_ev = resample(x_ev,0:0.01:x_ev.time(end));
    for i_x = 1:length(x_ev.Data)
        xy_c_t(i_x,:) = forward_kinematics(x_ev.Data(i_x,1),x_ev.Data(i_x,2),0,1);
    end
    hold on, plot(-xy_c_t(:,2),xy_c_t(:,1),'-', 'linewidth', 1, 'color', [0.40,0.40,0.40])
    
    hold off
    axis equal
    
    set(gca,'FontSize',18)
    xlabel('$x$','interpreter','latex','FontSize',22)
    ylabel('$y$','interpreter','latex','FontSize',22)
    set(gca,'TickLabelInterpreter','latex')
    grid on, box on

end