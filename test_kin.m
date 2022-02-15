%%%%%%%%%%%%%%%% Test Soft Kinematics %%%%%%%%%%%%%
clear all
close all
clc

%% Add functions path
addpath("my_functions");

%% Parameters
L = 1;
D = 0.05;

%% Curvature Trajectory
theta_step = 1;
% theta0 = -6*pi:theta_step:6*pi;
theta1 = -6*pi:theta_step:6*pi;

theta0 = zeros(1, length(theta1));
theta_r = (pi/6)*ones(1, length(theta0));

theta = [theta_r; theta0; theta1];

%% Simulate
figure
for i = 1:length(theta0)
   plot_Rsoft(theta(:, i), L, D, 0.1)
%    axis auto
   drawnow
end

