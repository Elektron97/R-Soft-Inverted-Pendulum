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
g = 9.81;
beta_r = 0.5;
beta = 0.1;
k = 1;
theta_init = [0; pi/4; -pi/4];
theta_dot_init = zeros(3, 1);

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
% load_system('master_thesis.slx');
% result=sim('master_thesis.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results

load_system('pick_place.slx');
result=sim('pick_place.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results

%% Plot and Animation
figure
plot(result.simout.time, result.simout.data)
grid on
xlabel("Time [s]");
ylabel("\theta [rad]");

% % Add desired trajectory

hold on
plot(result.simout2.time, result.simout2.data);
hold off
legend("\theta_r","\theta_0", "\theta_1", "\theta_{rd}");
title("Controlled R-SIP: Pick and Place");
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

% close all
% figure
% plot_Rsoft(result.simout.data(930, :), L, D, 0.1)
% plot_Rsoft([1.7; 3.1795; -3.81139], L, D, 0.1)
% hold on
% plot_Rsoft(result.simout.data(end, :), L, D, 0.1)
% hold off
%% Rec Video
rec = false;
close all

if rec
    v = VideoWriter("RSoft_Pendulum_PP.avi");
    open(v);
end

fig = figure;
fig.WindowState = 'fullscreen';
for i = 1:length(result.simout.time)
    subplot(2, 2, [1 3])
    plot_Rsoft([result.simout.data(i, 1); result.simout.data(i, 2); result.simout.data(i, 3)], L, D, 0.1)
   
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
