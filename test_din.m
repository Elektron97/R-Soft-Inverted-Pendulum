%%%%%%%%%%% Test Dinamica %%%%%%%%
clear all
close all
clc

addpath("my_functions");
addpath("Della Santina");
addpath("origin_soft_pendulum");
%% Parameters
L = 1;
D = 0.1;

m = 2;
g = 9.81;
% global k beta
beta = 0.1;
k = 1;
theta_init = [0; pi/4; -pi/4];
% theta_init = zeros(3, 1);
theta_dot_init = zeros(3, 1);

end_time = 30;
time_step = 0.01;
time = 0:time_step:end_time;

theta = theta_init;
theta_dot = theta_dot_init;

tau_r = 0;

%% Dynamics Matrix Validation
% %Perche' le mie vengono la meta'????????
% myB = originInertiaMatrix(theta(2), theta(3), 2*m, L, D);
% Bs = B_fcn(theta(2:3));
% Br = inertiaMatrix(theta(1), theta(2), theta(3), 2*m, L, D);
% 
% %Perche la gravita al contrarioooooooooooooo
% myG = originGravityVector(theta(2), theta(3), 2*m, -g, L, D);
% Gs = G_fcn(theta(2:3));
% Gr = gravityVector(theta(1), theta(2), theta(3), 2*m, -g, L, D);
% 
% myC = originCoriolisMatrix(theta(2), theta(3), theta_dot(2), theta_dot(3), 2*m, L, D);
% Cs = C_fcn(theta(2:3), theta_dot(2:3));
% Cr = coriolisMatrix(theta(1), theta(2), theta(3), theta_dot(1), theta_dot(2), theta_dot(3), 2*m, L, D);

figure
for i = 1:length(time)
    plot_Rsoft(theta, L, D, 0.1);
    drawnow
    
    %% Della Santina Simulation
%     ddx = f_fcn(theta(2:3), theta_dot(2:3), 0);
    
    %% My Simulation
    myB = originInertiaMatrix(theta(2), theta(3), m, L, D);
    myG = originGravityVector(theta(2), theta(3), m, -g, L, D);
    myC = originCoriolisMatrix(theta(2), theta(3), theta_dot(2), theta_dot(3), m, L, D);
    myD = originDampingMatrix(beta);
    myK = originElasticMatrix(k);
    
    ddx = -myB\(myC*theta_dot(2:3) + myG + myK*theta(2:3) + myD*theta_dot(2:3));
    
    theta_dot = theta_dot + time_step*[0;ddx];
    theta = theta + time_step*theta_dot;
end