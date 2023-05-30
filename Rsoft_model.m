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

% Consider elasticity in the revolute joint
syms k_r real

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here, Della Santina writes spatial coordinates, supposing the tangential%
% versor is: t(s) = [-sin(alpha); cos(alpha)].                            %
% In my notation, it has to be t(s) = [cos(alpha); sin(alpha)].           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % My version t(s) = [cos(alpha); sin(alpha)]
% x_s = L*(cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) ...
%      + sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));
% 
% y_s = -L*(sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) ...
%      - cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));

% Della Santina's version t(s) = [cos(alpha); sin(alpha)]
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

%% Forward Kinematics Function
alpha_I = theta_r + alpha;
pos = simplify(p_sd(1:2));

% matlabFunction([pos; alpha_I], 'File', 'forwardRSIP', 'Vars', [theta; L; D; d; s], 'Outputs', {'fwdKin'});

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

B =simplify( int( int(rho*(J_sd')*J_sd, d, [-0.5 0.5]), s, [0 1]) );
disp("Matrice di Inerzia Calcolata!");

%% Gravity Vector
gravity_field = int( int(rho*g*(sin(phi)*p_sd(1) + cos(phi)*p_sd(2)), d, [-0.5 0.5]), s, [0 1]);

G = sym(zeros(length(theta), 1));
for i = 1:length(theta)
    G(i) = simplify(diff(gravity_field, theta(i)));
end

disp("Gravità Calcolata!");

%% Elastic and Damping
H2 = henkelMatrix(2);

% including stiffness of revolute joint
K = blkdiag(k_r, k*H2);

Damp = blkdiag(beta_r, beta*H2);
%% Coriolis
C = christoffel(B, theta, theta_dot);
disp("Matrice di Coriolis Calcolata!");

%% Equilibria
potential = simplify(G + K*theta);

%% Equilibria: Varying Stiffness
% % To compute equilibria with numerical solutions
% step_k = 0.5;
% k_values = 1e-5:step_k:10;
% 
% % Number of trials
% n_try = 10;
% 
% % Autonomous System
% tau_r = 0;
% 
% for i = 1:length(k_values)
%     for j = 1:length(k_values)
%         equilibria_equation = simplify(subs(potential, [m; g; k_r; k; L; D], [1; 9.81; k_values(i); k_values(j); 1; 0.1])) == [1; 0; 0]*tau_r;
% 
%         % Numerical Solutions
%         for h = 1:n_try
%             solutions = vpasolve(equilibria_equation, theta, 'Random', true);
% 
%             if(isempty(solutions.theta_r))
%                 equilibria{i, j}(h, 1) = nan;
%                 equilibria{i, j}(h, 2) = nan;
%                 equilibria{i, j}(h, 3) = nan;
%             else
%                 equilibria{i, j}(h, 1) = double(solutions.theta_r);
%                 equilibria{i, j}(h, 2) = double(solutions.theta0);
%                 equilibria{i, j}(h, 3) = double(solutions.theta1);
%             end
% 
%             % Display to monitoring
%             disp("Iteration n " + num2str(h) + " | k = " + num2str(k_values(j)) + " | k_r = " + num2str(k_values(i)));
%         end
%     end
% end
% 
% save("equilibria_stiffness.mat", "equilibria");

%% Equilibria: Varying Torque
% % To compute equilibria with numerical solutions
% step_tau = 0.5;
% 
% tau_r = -10:step_tau:10;
% % rigid = 1e-5:0.5:10;
% n_try = 100;
% for i = 1:length(tau_r)
%     equilibria_equation = simplify(subs(potential, [m; g; k; L; D], [1; 9.81; 1; 1; 0.1])) == [1; 0; 0]*tau_r(i);
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
% 
% save("equilibria_tau.mat", "equilibria");

%% Equilibria: Zero Dynamics
% % Zero Dynamics (phi)
% % To compute equilibria with numerical solutions
% step_phi = pi/100;
% 
% phi = -pi:step_phi:pi;
% 
% n_try = 10;
% for i = 1:length(phi)
%     equilibria_equation = eval(simplify(subs(potential(2:3), [m; g; k; L; D; theta_r], [1; 9.81; 1; 1; 0.1; phi(i)])) == [0; 0]*0);
% 
% 
% 
%     for j = 1:n_try
%         solutions = vpasolve(equilibria_equation, theta, 'Random', true);
%         disp("Equilibrium for phi =" + num2str(phi(i)) + " | Tentative: " + num2str(j));
% 
%         if(isempty(solutions.theta0))
%             equilibria{i}(j, 1) = nan;
%             equilibria{i}(j, 2) = nan;
%         else
%             equilibria{i}(j, 1) = double(solutions.theta0);
%             equilibria{i}(j, 2) = double(solutions.theta1);
%         end
%     end
% end
% 
% % save("equilibria_phi.mat", "equilibria");

%% Computing Stiffness Matrix to Stability Analysis
for i = 1:length(potential)
    for j = 1: length(theta)
        Stiff_Mat(i, j) = diff(potential(i), theta(j));
    end
end

% matlabFunction(Stiff_Mat, 'File', 'stiffMatrix', 'Vars', [theta; m; g; k; L; D], 'Outputs', {'St_Mat'});
disp("Stiffness Matrix computed");

% K_equilibria = eval(eval(subs(Stiff_Mat, theta, [0; 0; 1e-9])))

%% State Space
% x1 = theta;
% x2 = theta_dot;
% 
% x = [x1; x2];
% inv_B = inv(B);
% S = [1 0 0]';
% 
% F = [x2; -inv_B*(C*x2 + G + K*x1 + D*x2)];
% G = [zeros(3, 1); inv_B*S];

%% Linearization
% Avoid recomputing inverse of inertia matrix
tic
disp("Computing Linearized System...");
inv_B = inv(B);

% A Matrix
A21 = -inv_B*(Stiff_Mat);
A22 = -inv_B*Damp;
A_lin = simplify(subs([zeros(3, 3), eye(3); A21, A22], theta_dot, zeros(3, 1)));
disp("Symbolic Linearization computed.");
toc

%% Save Functions
if(save_function)
    matlabFunction(B, 'File', 'inertiaMatrix', 'Vars', [theta; m; L; D], 'Outputs', {'B'});
    matlabFunction(G, 'File', 'gravityVector', 'Vars', [theta; m; g; L; D], 'Outputs', {'G'});
    matlabFunction(K, 'File', 'elasticMatrix', 'Vars', k, 'Outputs', {'K'});
    matlabFunction(Damp, 'File', 'dampingMatrix', 'Vars', [beta; beta_r], 'Outputs', {'D'});
    matlabFunction(C, 'File', 'coriolisMatrix', 'Vars', [theta; theta_dot; m; L; D], 'Outputs', {'C'});
    matlabFunction(Stiff_Mat, 'File', 'stiffMatrix', 'Vars', [theta; m; g; k; L; D], 'Outputs', {'St_Mat'});
    
% % Doesn't work, thanks matlab
%     matlabFunction(B, 'File', 'uniformInertiaMatrix', 'Vars', [theta; m; L; D], 'Outputs', {'B'});
%     matlabFunction(G, 'File', 'uniformGravityVector', 'Vars', [theta; m; g; L; D], 'Outputs', {'G'});
%     matlabFunction(K, 'File', 'uniformElasticMatrix', 'Vars', k, 'Outputs', {'K'});
%     matlabFunction(Damp, 'File', 'uniformDampingMatrix', 'Vars', [beta; beta_r], 'Outputs', {'D'});
%     matlabFunction(C, 'File', 'uniformCoriolisMatrix', 'Vars', [theta; theta_dot; m; L; D], 'Outputs', {'C'});
end
