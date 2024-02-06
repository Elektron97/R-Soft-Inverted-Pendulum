%%%%%%%%%%% Final Simulations %%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");
addpath("Della Santina");
addpath("origin_soft_pendulum");
%% Parameters
L = 1;
D = 0.1;

m = 1;
m_obj = 0.0;
g = 9.81;
beta_r = 0.5;
beta = 0.1;
k = 1;
theta_init = [0; pi/4; -pi/4];
theta_dot_init = zeros(3, 1);

% Useful colors
gray_col = [0.7, 0.7, 0.7];
htmlGray = [128 128 128]/255;
purple_col = [0.4940 0.1840 0.5560];
orange_col = [0.8500 0.3250 0.0980];
green_col = [0.4660 0.6740 0.1880];
celeste_figo = [0 0.4470 0.7410];
cool_yellow = [0.9290 0.6940 0.1250];
cool_red = [0.6350 0.0780 0.1840];
black = [0 0 0];

%% LQR: optimal gains PD
A_lin = [0 1; 0 0];
B_lin = [0; -1];
Q = 100*eye(2);
R = 100;

% lqr() gives K | A - BK
[K, ~, ~] = lqr(A_lin, B_lin, Q, R);
Kp = -K(1);
Kd = -K(2);

%% Simulation
is_pp = false;
if is_pp
    load_system('pick_place.slx');
    result=sim('pick_place.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results
else
    load_system('master_thesis.slx');
    result=sim('master_thesis.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results
end

%% Plot and Animation
figure
plot(result.simout.time, result.simout.data, 'LineWidth', 2.0)
grid on
xlabel("Time [s]");
ylabel("\theta [rad]");

% % Add desired trajectory

hold on
plot(result.simout2.time, result.simout2.data, '--', 'LineWidth', 2.0);
hold off
legend("\theta_r","\theta_0", "\theta_1", "\theta_{rd}");
% title("Controlled R-SIP: Pick and Place");
ylim([-15, 15])

% % Phase Space

figure
plot3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1))
title("Space of Phase")
grid on
xlabel("\theta_0")
ylabel("\theta_1")
zlabel("\theta_r")

% % Add initial and final position

hold on
plot3(result.simout.data(1, 2), result.simout.data(1, 3), result.simout.data(1, 1), 'rx')
plot3(result.simout.data(end, 2), result.simout.data(end, 3), result.simout.data(end, 1), 'gx')
hold off


% Add quiver3
hold on
quiver3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1), ...
        result.simout3.data(:, 2), result.simout3.data(:, 3), result.simout3.data(:, 1), 'Color', [0.9290 0.6940 0.1250])
hold off
axis equal

% % Actuation
% figure
% plot(result.simout1.time, result.simout1.data)
% grid on
% xlabel("Time [s]");
% ylabel("\tau [N m]");
% title("Actuation");

%% Multiple Frame Plot
% close all
% 
% % Insert skip of robot frame configuration
% skip_frame_robot = 5;
% % skip_frame_quiver = 7;
% 
% if is_pp
%     phase_time = [1, 739, 928, length(result.simout.time)];
%     phase_color = {green_col, orange_col, cool_yellow, celeste_figo};
% 
%     %%%%% Cicle on Phases
%     for phase=1:3
% 
%         figure
%         hold on
%         % Initial Position
%         plot_Rsoft(result.simout.data(phase_time(phase), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase})
%         % Final Position
%         plot_Rsoft(result.simout.data(phase_time(phase+1), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase+1})
%         for i = phase_time(phase):phase_time(phase+1)
%             % skip_frame
%             % rewrite (i % skip_frame_robot == 0) in matlab
%             if(mod(i, skip_frame_robot) == 0)
%                 plot_Rsoft(result.simout.data(i, :), L, D, plot_frame=false, plot_thick=false, color=gray_col, backbone_thick=1)
%             end
%             % Trajectory always continous
%             tip_traj(:, i) = fwdKinRSIP(result.simout.data(i, :), 1, 0, L, D);
%             % tip_vel(:, i) = jacobianSoft(result.simout.data(i, 1), ...
%             %                                 result.simout.data(i, 2), ...
%             %                                 result.simout.data(i, 3), ...
%             %                                 1, 0, L, D)*result.simout3.data(i, :)';
%         end
%         % Initial Position
%         plot_Rsoft(result.simout.data(phase_time(phase), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase})
%         % Final Position
%         plot_Rsoft(result.simout.data(phase_time(phase+1), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase+1})
%         hold on
%         % Tip Trajectory
%         plot3(tip_traj(1, phase_time(phase):phase_time(phase+1)), ...
%             tip_traj(2, phase_time(phase):phase_time(phase+1)), ...
%             zeros(1, length(phase_time(phase):phase_time(phase+1))), ...
%             'LineWidth', 2, 'Color', htmlGray)
%         %Tip Velocity: all points
%         % quiver3(tip_traj(1, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
%         %         tip_traj(2, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
%         %         zeros(1, length(phase_time(phase):skip_frame_quiver:phase_time(phase+1))), ...
%         %         tip_vel(1, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
%         %         tip_vel(2, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
%         %         zeros(1, length(phase_time(phase):skip_frame_quiver:phase_time(phase+1))), ...
%         %         'Color', '#25283D', 'LineWidth', 1.0)
%         hold off
%         xlabel("x [m]")
%         ylabel("y [m]")
%         % title("Pick and Place: Phase " + num2str(phase))
%         legend("Initial Config. of Phase " + num2str(phase), "Final Config. of Phase " + num2str(phase))
% 
%         clear tip_traj
%         clear tip_vel
%     end
% else
%     figure
%     hold on
%     for i =1:length(result.simout.time)
%         plot_Rsoft(result.simout.data(i, :), L, D, plot_frame=false, plot_thick=false, color=gray_col, backbone_thick=1)
%         tip_traj(:, i) = fwdKinRSIP(result.simout.data(i, :), 1, 0, L, D);
%     end
%     % Initial Position
%     plot_Rsoft(result.simout.data(1, :), L, D, plot_frame=false, plot_thick=false, color=green_col)
%     % Final Position
%     plot_Rsoft(result.simout.data(end, :), L, D, plot_frame=false, plot_thick=false, color=orange_col)
%     hold on
%     % Tip Trajectory
%     plot3(tip_traj(1, :), tip_traj(2, :), zeros(1, length(result.simout.time)), ...
%         'LineWidth', 2, 'Color', htmlGray)
%     hold off
%     xlabel("x [m]")
%     ylabel("y [m]")
% end
%% Rec Video
% Tip 
rec = false;
close all

if rec
    v = VideoWriter("RSoft_Pendulum_PP.avi");
    open(v);
end

fig = figure;
% fig.WindowState = 'fullscreen';
for i = 1:length(result.simout.time)
    subplot(2, 2, [1 3])
        cla reset % avoid multiple frame
        plot_Rsoft([result.simout.data(i, 1); result.simout.data(i, 2); result.simout.data(i, 3)], L, D);
        hold off
   
    subplot(2, 2, 2)
        plot(result.simout.time(1:i), result.simout.data(1:i, :))
        hold on
         plot(result.simout2.time(1:i), result.simout2.data(1:i, :))
        hold off
        grid on
        legend("\theta_r","\theta_0", "\theta_1", "\theta_rd");
        xlabel("Time [s]");
        ylabel("\theta [rad]");
        title("Trajectory");

    
    subplot(2, 2, 4)
    %     plot(result.simout1.time(1:i), result.simout1.data(1:i))
    %     grid on
    %     xlabel("Time [s]");
    %     ylabel("\tau [N m]");
    %     title("Actuation");

        plot3(result.simout.data(1:i, 2), result.simout.data(1:i, 3), result.simout.data(1:i, 1))
        title("Phase Space")
        grid on
        xlabel("\theta_0")
        ylabel("\theta_1")
        zlabel("\theta_r")
    
    % % Add initial and final position
    
    %     hold on
    %     plot3(result.simout.data(1, 2), result.simout.data(1, 3), result.simout.data(1, 1), 'rx')
    %     plot3(result.simout.data(end, 2), result.simout.data(end, 3), result.simout.data(end, 1), 'gx')
    %     hold off
    
    
    % Add quiver3
        hold on
        quiver3(result.simout.data(1:i, 2), result.simout.data(1:i, 3), result.simout.data(1:i, 1), ...
                result.simout3.data(1:i, 2), result.simout3.data(1:i, 3), result.simout3.data(1:i, 1))
        hold off
        axis equal

    if rec
        frame = getframe(gcf);
        writeVideo(v, frame);
    else
        drawnow
    end

end

if rec
    close(v);
end
