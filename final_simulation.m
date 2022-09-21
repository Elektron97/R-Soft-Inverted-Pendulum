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
load_system('master_thesis.slx');
result=sim('master_thesis.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results

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
title("Controlled R-SIP: Step traj.");
ylim([-15, 15])

figure
plot3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1))
title("Spazio delle fasi")
grid on
xlabel("\theta_0")
ylabel("\theta_1")
zlabel("\theta_r")

% figure
% plot(result.simout1.time, result.simout1.data)
% grid on
% xlabel("Time [s]");
% ylabel("\tau [N m]");
% title("Actuation");

%% Add equilibria in phase space
zero_equilibria
for i = 1:length(phi)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)     
        if(stability{i}(j))
            plot(phi(i), filtered_equilibria{i}(j, 1), 'bx')
        else
            plot(phi(i), filtered_equilibria{i}(j, 1), 'rx')
        end 
        
    end
end

%% Rec Video
% v = VideoWriter("RSoft_Pendulum_step");
% open(v);
% 
% fig = figure;
% fig.WindowState = 'fullscreen';
% for i = 1:length(result.simout.time)
%     subplot(2, 2, [1 3])
% %     plot_Rsoft([0; result.simout.data(i, 1); result.simout.data(i, 2)], L, D, 0.1)
%     plot_Rsoft([result.simout.data(i, 1); result.simout.data(i, 2); result.simout.data(i, 3)], L, D, 0.1)
%    
%     subplot(2, 2, 2)
%     plot(result.simout.time(1:i), result.simout.data(1:i, :))
%     grid on
% %     legend("\theta_0", "\theta_1");
%     legend("\theta_r","\theta_0", "\theta_1");
%     xlabel("Time [s]");
%     ylabel("\theta [rad]");
%     title("Snap Trajectory");
%     
%     subplot(2, 2, 4)
%     plot(result.simout1.time(1:i), result.simout1.data(1:i))
%     grid on
%     xlabel("Time [s]");
%     ylabel("\tau [N m]");
%     title("Actuation");
% 
% %     drawnow
%     frame = getframe(gcf);
%     writeVideo(v, frame);
% end
% 
% close(v);