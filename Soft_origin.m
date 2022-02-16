%%%%%%%% Soft Inverted %%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");

save_functions = false;
%% Declare Symbolic Variables
syms theta0 theta1 real
syms theta0_dot theta1_dot real
syms s d real
syms L D real

theta = [theta0; theta1];
theta_dot = [theta0_dot; theta1_dot];

syms m k g beta real

phi = 0;
%% Affine Curvature
K = theta0 + theta1*s;

%% Orientation of SoR {Ss} w.r.t. {S0}
alpha = int(K, s, 0, s);

fresn_sin1 = fresnels( (theta0 + s*theta1) *  sqrt(1/(pi*theta1)) );
fresn_sin2 = fresnels(theta0 *  sqrt(1/(pi*theta1)) );

fresn_cos1 = fresnelc((theta0 + s*theta1) * sqrt(1/(pi*theta1)));
fresn_cos2 = fresnelc(theta0 * sqrt(1/(pi*theta1)) );

x_s = L*(sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) ...
     - cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));
 
y_s = L*(cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) ...
     + sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));
 
p_s = simplify([x_s; y_s; 0]);

T0s = [my_rot(alpha, 'z') p_s; zeros(1, 3) 1] ;

%% Point on Thickness D
p_sd_hom = simplify(T0s*[d*D; 0; 0; 1]);
p_sd = p_sd_hom(1:3);

%% Differential Kinematics
J_sd = sym(zeros(2, length(theta)));

%i-joint variables
for i = 1:length(theta)
        J_sd(1, i) = diff(p_sd(1), theta(i));
        J_sd(2, i) = diff(p_sd(2), theta(i));
end

J_sd = simplify(J_sd);

twist_s = J_sd * theta_dot;

%% Inertia Matrix
rho = m*dirac(s-1);
B = simplify(int( int(rho*(J_sd')*J_sd, d, [-0.5 0.5]), s, [0 1]));

%% Gravity Vector
gravity_field = int( int(rho*g*(sin(phi)*p_sd(1) + cos(phi)*p_sd(2)), d, [-0.5 0.5]), s, [0 1]);

G = sym(zeros(length(theta), 1));
for i = 1:length(theta)
    G(i) = simplify(diff(gravity_field, theta(i)));
end

%% Elastic and Damping
H2 = henkelMatrix(2);

K = k*H2;
Damp = beta*H2;
%% Coriolis
C = christoffel(B, theta, theta_dot);

%% Equilibri
%Stiffness Matrix
potential = G + K*theta;

for i = 1:length(potential)
    for j = 1: length(theta)
        Stiff_Mat(i, j) = diff(potential(i), theta(j));
    end
end

Stiff_Mat_eq = eval(eval(subs(Stiff_Mat, theta, (1e-5)*ones(length(theta), 1))));

% %Definita positiva
% det(eval(eval(Stiff_Mat_eq(1, 1)))) > 0
% det(eval(eval(Stiff_Mat_eq(1:2, 1:2)))) > 0
% det(eval(eval(Stiff_Mat_eq))) > 0

%% Save Functions
if(save_functions)
    matlabFunction(B, 'File', 'originInertiaMatrix', 'Vars', [theta; m; L; D], 'Outputs', {'B'});
    matlabFunction(G, 'File', 'originGravityVector', 'Vars', [theta; m; g; L; D], 'Outputs', {'G'});
    matlabFunction(K, 'File', 'originElasticMatrix', 'Vars', k, 'Outputs', {'K'});
    matlabFunction(Damp, 'File', 'originDampingMatrix', 'Vars', beta, 'Outputs', {'D'});
    matlabFunction(C, 'File', 'originCoriolisMatrix', 'Vars', [theta; theta_dot; m; L; D], 'Outputs', {'C'});
end

