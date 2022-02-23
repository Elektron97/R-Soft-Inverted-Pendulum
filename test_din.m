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

m = 1;
g = 9.81;
% global k beta
beta_r = 0.5;
beta = 0.1;
k = 4;
theta_init = [0; pi/4; -pi/4];
% theta_init = zeros(3, 1);
theta_dot_init = zeros(3, 1);

end_time = 30;
time_step = 0.01;
time = 0:time_step:end_time;

theta = theta_init;
theta_dot = theta_dot_init;

tau_r = 0;

%% Simulation

figure
for i = 1:length(time)
    plot_Rsoft(theta, L, D, 0.1);
    drawnow
    
    %% Della Santina Simulation
%     ddx = f_fcn(theta(2:3), theta_dot(2:3), 0);
%     theta_dot = theta_dot + time_step*[0;ddx];
%     theta = theta + time_step*theta_dot;
    
    %% My Simulation    
%     %Origin Soft Inverted Pendulum
    ddx = Soft_dynamics(theta(2:3), theta_dot(2:3), 0, m, g, L, D, k, beta);
    theta_dot = theta_dot +time_step*[0; ddx];
    theta = theta + time_step*theta_dot;
    
    %R-Soft Inverted Pendulum
%     ddx = RSoft_dynamics(theta, theta_dot, tau_r, m, g, L, D, k, beta);
%     
%     theta_dot = theta_dot + time_step*ddx;
%     theta = theta + time_step*theta_dot;
end