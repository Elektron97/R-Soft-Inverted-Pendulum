%%%%%%%%%%% R-Soft Inverted Pendulum %%%%%%%%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");

save_function = false;
%% Declare Symbolic Variables
syms theta_r theta0 theta1 real
syms theta_r_dot theta0_dot theta1_dot real
syms s d real
syms L D real

theta = [theta_r; theta0; theta1];
theta_dot = [theta_r_dot; theta0_dot; theta1_dot];

syms m k g beta beta_r real

phi = 0;
%% Orientation of SoR {S0} w.r.t. {I}
Ri0 = my_rot(theta_r, 'z');
Ti0 = blkdiag(Ri0, 1);

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

J_sd = simplify(J_sd);

twist_s = J_sd * theta_dot;
%% Inertia Matrix
%Fix dirac bug
rho = 2*m*dirac(s-1);
% rho = m;

B = simplify(int( int(rho*(J_sd')*J_sd, d, [-0.5 0.5]), s, [0 1]));

%% Test B rotazionale
% A = [s; s^2/2];
% A_r = [1; A];
% 
% B_rot = simplify(int( int(rho*(p_sd(1)^2 + p_sd(2)^2)*A_r*(A_r'), d, [-0.5 0.5]), s, [0 1]));
% % matlabFunction(B_rot, 'File', 'rotInertiaMatrix', 'Vars', [theta; m; L; D], 'Outputs', {'B'});
%% Gravity Vector
gravity_field = int( int(rho*g*(sin(phi)*p_sd(1) + cos(phi)*p_sd(2)), d, [-0.5 0.5]), s, [0 1]);

G = sym(zeros(length(theta), 1));
for i = 1:length(theta)
    G(i) = simplify(diff(gravity_field, theta(i)));
end

%% Elastic and Damping
H2 = henkelMatrix(2);

K = k*blkdiag(0, H2);
Damp = blkdiag(beta_r, beta*H2);
%% Coriolis
C = christoffel(B, theta, theta_dot);

%% Equilibria
% potential = simplify(G + K*theta);
% 
% step_tau = 0.1;
% 
% % tau_r = -5:step_tau:5;
% rigid = 1e-5:0.5:10;
% tau_r = 0;
% n_try = 100;
% for i = 1:length(rigid)
%     equilibria_equation = simplify(subs(potential, [m; g; k; L; D], [1; 9.81; rigid(i); 1; 0.1])) == [1; 0; 0]*tau_r;
% 
% 
% 
%     for j = 1:n_try
%         solutions = vpasolve(equilibria_equation, theta, 'Random', true);
% 
%         if(isempty(solutions.theta_r))
%             equilibria{i}(j, 1) = nan;
%             equilibria{i}(j, 2) = nan;
%             equilibria{i}(j, 3) = nan;
%         else
%             equilibria{i}(j, 1) = double(solutions.theta_r);
%             equilibria{i}(j, 2) = double(solutions.theta0);
%             equilibria{i}(j, 3) = double(solutions.theta1);
%         end
%     end
% end

% disp("Equilibria: theta_R = " + num2str(wrapToPi(double(equilibria.theta_r))) + ...
%         " theta0 = " + num2str(wrapToPi(double(equilibria.theta0))) + " theta1 = " + ...
%             num2str(wrapToPi(double(equilibria.theta1))));

% for i = 1:length(potential)
%     for j = 1: length(theta)
%         Stiff_Mat(i, j) = diff(potential(i), theta(j));
%     end
% end
% 
% Stiff_Mat_eq1 = eval(eval(subs(simplify(Stiff_Mat), theta, 1e-5*ones(3, 1))));
% Stiff_Mat_eq2 = eval(eval(subs(simplify(Stiff_Mat), theta, [pi; 1e-5; 1e-5])));
%% Save Functions
if(save_function)
    matlabFunction(B, 'File', 'inertiaMatrix', 'Vars', [theta; m; L; D], 'Outputs', {'B'});
    matlabFunction(G, 'File', 'gravityVector', 'Vars', [theta; m; g; L; D], 'Outputs', {'G'});
    matlabFunction(K, 'File', 'elasticMatrix', 'Vars', k, 'Outputs', {'K'});
    matlabFunction(Damp, 'File', 'dampingMatrix', 'Vars', [beta; beta_r], 'Outputs', {'D'});
    matlabFunction(C, 'File', 'coriolisMatrix', 'Vars', [theta; theta_dot; m; L; D], 'Outputs', {'C'});
end
