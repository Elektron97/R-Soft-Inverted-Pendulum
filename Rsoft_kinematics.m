%%%%%%%%%%% R-Soft Inverted Pendulum %%%%%%%%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");

%% Declare Symbolic Variables
syms theta_r theta0 theta1 real
syms theta_r_dot theta0_dot theta1_dot real
syms s d real
syms L D real

theta = [theta_r; theta0; theta1];
theta_dot = [theta_r_dot; theta0_dot; theta1_dot];
%% Orientation of SoR {S0} w.r.t. {I}
Ri0 = my_rot(theta_r, 'z');
Ti0 = blkdiag(Ri0, 1);

%% Affine Curvature
K = theta0 + theta1*s;

%% Orientation of SoR {Ss} w.r.t. {S0}
alpha = int(K, s, 0, s);

fresn_sin1 = fresnels((theta0 + s*theta1)/(sqrt(pi*theta1)));
fresn_sin2 = fresnels(theta0/sqrt(pi*theta1));

fresn_cos1 = fresnelc((theta0 + s*theta1)/(sqrt(pi*theta1)));
fresn_cos2 = fresnelc(theta0/sqrt(pi*theta1));

x_s = L*(sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) ...
     - cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));
 
y_s = L*(cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) ...
     + sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));
 
p_s = simplify([x_s; y_s; 0]);

T0s = [my_rot(alpha, 'z') p_s; zeros(1, 3) 1] ;

%% Point on Thickness D
p_sd_hom = simplify(Ti0*T0s*[d*D; 0; 0; 1]);
p_sd = p_sd_hom(1:3);

%Validation
limit(p_sd, theta1, 0)

%% Differential Kinematics
J_sd = sym(zeros(2, length(theta)));

%i-joint variables
for i = 1:length(theta)
    %j: 2D coordinates
    for j=1:2
        J_sd(j, i) = diff(p_sd(j), theta(i));
    end
end

twist_s = J_sd * theta_dot;



