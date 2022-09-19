%%%%%%%%%%%%% Test filtrazione %%%%%%%%
clear all
close all
clc

addpath("my_functions");
%% Models
% syms x1 x2 real
% 
% f = [x2^2; 0];
% g = [0 1]';
% x = [x1; x2];
% 
% filtration([f g], g, x)
% 
% syms x1 x2 x3 real
% 
%  f = [2*x2 1; 1 0; 0 x2];
%  g = [];
% 
%  x = [x1 x2 x3]';
% 
%  filtration([f g], g, x)

% syms x y theta real
% 
% f = [];
% g = [cos(theta) 0; sin(theta) 0; 0 1];
% dh = [x y 0];
% 
% filtration([f g], g, [x y theta]')
% rowFiltration([f g], dh, [x y theta]')

% syms x1 x2 x3 real
% x = [x1 x2 x3];
% 
% f = [x1*x2 0 0]';
% g = [x3 0 0]';
% h = x1;
% 
% rowFiltration([f g], jacobian(h, x), x)

% syms x y theta real
% syms d real
% state = [x y theta]';
% f = [];
% g = [cos(theta) 0; sin(theta) 0; 0 1];
% h1 = pi - theta + atan2(y, x);
% h2 = pi - theta - atan2(y - d, x);
% 
% dh = [jacobian(h1, state); jacobian(h2, state)];
% 
% % filtration([f g], g, [x y theta]')
% rowFiltration([f g], dh, [x y theta]')

syms x y phi x_dot y_dot phi_dot m g L J real

f = [x_dot; y_dot; phi_dot; 0; -g; 0];
g1 = [0; 0; 0; -sin(phi)/m; cos(phi)/m; -L/J];
g2 = [g1(1:end-1); -g1(end)];
g = [g1 g2];

filtration([f g], g, [x y phi x_dot y_dot phi_dot]')


