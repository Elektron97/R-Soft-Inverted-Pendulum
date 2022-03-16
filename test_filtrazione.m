%%%%%%%%%%%%% Test filtrazione %%%%%%%%
clear all
close all
clc

addpath("my_functions");
%% Symbolic variables
% syms x1 x2 real
% 
% f = [x2^2; 0];
% g = [0 1]';
% x = [x1; x2];
% 
% filtration([f g], g, x)

% syms x1 x2 x3 real
% 
%  f = [2*x2 1; 1 0; 0 x2];
%  g = [];
% 
%  x = [x1 x2 x3]';
% 
%  filtration([f g], g, x)

syms x y theta real

f = [];
g = [cos(theta) 0; sin(theta) 0; 0 1];

filtration([f g], g, [x y theta]')