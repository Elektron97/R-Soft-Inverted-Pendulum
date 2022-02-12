%%%%%%%%%%% Test Dinamica %%%%%%%%
clear all
close all
clc

addpath("my_functions");

%% Parameters
L = 1;
D = 0.1;

m = 1;
g = 9.81;
beta = 0.1;
k = 1;
% theta_init = [0; pi/4; -pi/4];
theta_init = zeros(3, 1);
theta_dot_init = zeros(3, 1);

end_time = 0.1;
time_step = 0.001;
time = 0:time_step:end_time;

theta = theta_init;
theta_dot = theta_dot_init;

tau_r = 0;

% B = inertiaMatrix(theta(1), theta(2), theta(3), m, L, D);
% G = gravityVector(theta(1), theta(2), theta(3), m, g, L, D);
% C = coriolisMatrix(theta(1), theta(2), theta(3), theta_dot(1), theta_dot(2), theta_dot(3), m, L, D);
% K = elasticMatrix(k);
% D = dampingMatrix(beta);

figure
for i = 1:length(time)
    plot_Rsoft(theta, L, D, 0.1);
    drawnow
    
    %State Model
    x1 = theta;
    x2 = theta_dot;
    
    if abs(x1(2)) < 1e-5
        x1(2) = 1e-5;
    end

    if abs(x1(3)) < 1e-5
        x1(3) = 1e-5;
    end

    B = inertiaMatrix(theta(1), theta(2), theta(3), m, L, D);
    G = gravityVector(theta(1), theta(2), theta(3), m, g, L, D);
    C = coriolisMatrix(theta(1), theta(2), theta(3), theta_dot(1), theta_dot(2), theta_dot(3), m, L, D);
    K = elasticMatrix(k);
    Damp = dampingMatrix(beta);
    
    F = [x2; -B\(C*x2 + G + K*x1 + Damp*x2)];
    G = [zeros(3, 1); B\[1; 0; 0]];
    
    x_dot = F + G*tau_r;
    
    %Update
    theta = theta + time_step*x_dot(1:3);
    theta_dot = theta_dot + time_step*x_dot(4:6);
end



% %% State Model Function
% function x_dot = state_model(theta, theta_dot, tau_r, m, L, D, g, k, beta)
% x1 = theta;
% x2 = theta_dot;
% 
% B = inertiaMatrix(theta(1), theta(2), theta(3), m, L, D);
% G = gravityVector(theta(1), theta(2), theta(3), m, g, L, D);
% C = coriolisMatrix(theta(1), theta(2), theta(3), theta_dot(1), theta_dot(2), theta_dot(3), m, L, D);
% K = elasticMatrix(k);
% D = dampingMatrix(beta);
% 
% f = [x2; -B\(C*x2 + G + K*x1 + D*x2)];
% g = [zeros(3, 1); B\[1; 0; 0]];
% 
% x_dot = f + g*tau_r;
% end

