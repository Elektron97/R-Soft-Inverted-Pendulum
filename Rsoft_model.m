%%%%%%%%%%% R-Soft Inverted Pendulum %%%%%%%%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");

save_function = true;
%% Declare Symbolic Variables
syms theta_r theta0 theta1 real
syms theta_r_dot theta0_dot theta1_dot real
syms s d real
syms L D real

theta = [theta_r; theta0; theta1];
theta_dot = [theta_r_dot; theta0_dot; theta1_dot];

syms m k g beta real

phi = 0;
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

Tis = Ti0*T0s;

%% Point on Thickness D
p_sd_hom = simplify(Tis*[d*D; 0; 0; 1]);
p_sd = p_sd_hom(1:3);

%% Differential Kinematics
J_sd = sym(zeros(2, length(theta)));

%i-joint variables
for i = 1:length(theta)
        J_sd(1, i) = diff(p_sd(1), theta(i));
        J_sd(2, i) = diff(p_sd(2), theta(i));
end

twist_s = J_sd * theta_dot;

%% Inertia Matrix
rho = m*dirac(s-1);
% rho = m;
B = int( int(rho*(J_sd')*J_sd, d, [-0.5 0.5]), s, [0 1]);

%% Gravity Vector
gravity_field = int( int(rho*g*(sin(phi)*p_sd(1) + cos(phi)*p_sd(2)), d, [-0.5 0.5]), s, [0 1]);

G = sym(zeros(length(theta), 1));
for i = 1:length(theta)
    G(i) = diff(gravity_field, theta(i));
end

%% Elastic and Damping
H2 = henkelMatrix(2);

K = k*blkdiag(0, H2);
Damp = beta*blkdiag(0, H2);
%% Coriolis
C = christoffel(B, theta, theta_dot);

%% Save Functions
if(save_function)
    matlabFunction(B, 'File', 'inertiaMatrix', 'Vars', [theta; m; L; D], 'Outputs', {'B'});
    matlabFunction(G, 'File', 'gravityVector', 'Vars', [theta; m; g; L; D], 'Outputs', {'G'});
    matlabFunction(K, 'File', 'elasticMatrix', 'Vars', k, 'Outputs', {'K'});
    matlabFunction(Damp, 'File', 'dampingMatrix', 'Vars', beta, 'Outputs', {'D'});
    matlabFunction(C, 'File', 'coriolisMatrix', 'Vars', [theta; theta_dot; m; L; D], 'Outputs', {'C'});
end
