global L k x_0 dx_0 beta tau

global k beta
global m g L

addpath('automatically_generated')

m = 1;
g = 9.81;
k = 0.1;
L = 1;
D = 0.1;
beta = .1;
tau = 0;

k/(m*g*L) > (13+2*sqrt(31))/60
k/(m*g*L) > 1/10

%  theta_0_0 = pi/4;
%  theta_0_1 = 0;
% dtheta_0_0 = 0;
% dtheta_0_1 = 0;

x_0 = [ pi/4; -pi/4 ];
% x_0 = [pi; eps];
% x_0 = [6.1881;   -6.8723];
% x_0 = [10.2863;  -15.2886];
dx_0 = [    0; 0 ];