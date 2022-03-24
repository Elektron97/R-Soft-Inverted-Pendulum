%%%%%%%%%%% Test Dinamica %%%%%%%%
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
% D = 0;

m = 1;
g = 9.81;
beta_r = 0.5;
beta = 0.1;
k = 1;
theta_init = [0; pi/4; -pi/4];
theta_dot_init = zeros(3, 1);

Kp = 5;
Kd = 2;

% thetaR_des = pi/4;
alpha_des = 0;

%% Adaptive initial Parameters
% pi_real = [m*L^2; m*g*L; k; beta];
% pi_init = [0.8; 9.7; 0.91; 0.15];
pi_init = zeros(4, 1);
%% Simulation
% load_system('R_soft_sim.slx');
% result=sim('R_soft_sim.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results

load_system('R_soft_sim2020b.slx');
result=sim('R_soft_sim2020b.slx', 'ReturnWorkspaceOutputs','on'); %simulate and extract results

%% Plot and Animation
figure
subplot(2, 1, 1)
plot(result.simout.time, result.simout.data)
grid on
% legend("\theta_0", "\theta_1");
legend("\theta_r","\theta_0", "\theta_1");
xlabel("Time [s]");
ylabel("\theta [rad]");

subplot(2, 1, 2)
plot(result.simout1.time, result.simout1.data)
grid on
xlabel("Time [s]");
ylabel("\tau [N m]");
title("Actuation");

% v = VideoWriter("RSoft_Pendulum_alpha");
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