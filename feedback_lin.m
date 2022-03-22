%%%%%%%%%%%% Feedback Linearization %%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Model
Rsoft_model

%% Redefine state variables
x1 = theta_r;
x2 = [theta0; theta1];
x3 = theta_r_dot;
x4 = [theta0_dot; theta1_dot];
x = [x1; x2; x3; x4];
%% Alpha Output
h = [1 1 1/2 0 0 0]*x;
LfH = jacobian(h, x)*F;
LgH = jacobian(h, x)*G;
Lf2H = jacobian(LfH, x)*F;
LgLfH = jacobian(LfH, x)*G;

Phi = [1 1 1/2 0 0 0; 0 0 0 1 1 1/2; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]*x;
diffPhi = jacobian(Phi, x);
% Verify that is a consistent change of variables
% det(diffPhi)

syms csi1 csi2 eta1 eta2 eta3 eta4 real
new_state = [csi1; csi2; eta1; eta2; eta3; eta4];

syms v real
u = -(Lf2H/LgLfH) + (1/LgLfH)*v;

%% Zero Dynamics
disp("Zero Dynamics")
f_zero = simplify(subs([eta3; eta4; subs([zeros(2, 4) eye(2, 2)]*F, x, inv(diffPhi)*new_state)], [csi1; csi2], zeros(2, 1)));
disp("Equilibria")
equilibria_equation = f_zero == zeros(4, 1);

n_try = 100;
for j = 1:n_try
    solutions = vpasolve(equilibria_equation, [eta1; eta2; eta3; eta4], 'Random', true);
    if(isempty(solutions.eta1))
        equilibria{i}(j, 1) = nan;
        equilibria{i}(j, 2) = nan;
        equilibria{i}(j, 3) = nan;
        equilibria{i}(j, 4) = nan;
    else
        equilibria{i}(j, 1) = double(solutions.eta1);
        equilibria{i}(j, 2) = double(solutions.eta2);
        equilibria{i}(j, 3) = double(solutions.eta3);
        equilibria{i}(j, 4) = double(solutions.eta4);
    end
end

save("equilibriaZeroAlpha.mat", "equilibria");

%% Cartesian Output
% % h = [1 1 1/2 0 0 0]*x;
% h = simplify(cos(theta_r)*subs(x_s, s, 1) - sin(theta_r)*subs(y_s, s, 1));
% LfH = jacobian(h, x)*F;
% LgH = jacobian(h, x)*G;
% Lf2H = jacobian(LfH, x)*F;
% LgLfH = jacobian(LfH, x)*G;
% 
% % Phi = [1 1 1/2 0 0 0; 0 0 0 1 1 1/2; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]*x;
% % diffPhi = jacobian(Phi, x);
% % % Verify that is a consistent change of variables
% % % det(diffPhi)
% % 
% % new_state = diffPhi*x;
% % csi1 = new_state(1);
% % csi2 = new_state(2);
% % eta1 = new_state(3);
% % eta2 = new_state(4);
% % eta3 = new_state(5);
% % eta4 = new_state(6);
% 
% syms v real
% u = -(Lf2H/LgLfH) + (1/LgLfH)*v;
% disp("Save control in a matlab Function");
% matlabFunction(u, 'File', 'xsFL', 'Vars', [theta; theta_dot; m; g; k; L; D; beta_r; beta; v], 'Outputs', {'xsFL'})
