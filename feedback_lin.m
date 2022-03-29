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
% h = [1 1 1/2 0 0 0]*x;
% LfH = jacobian(h, x)*F;
% LgH = jacobian(h, x)*G;
% Lf2H = jacobian(LfH, x)*F;
% LgLfH = jacobian(LfH, x)*G;
% disp("Secondo ordine calcolato");
% 
% % Lf3H = jacobian(Lf2H, x)*F;
% % Lg2LfH = jacobian(Lf2H, x)*G;
% 
% Phi = [1 1 1/2 0 0 0; 0 0 0 1 1 1/2; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]*x;
% diffPhi = jacobian(Phi, x);
% % Verify that is a consistent change of variables
% % det(diffPhi)
% 
% syms csi1 csi2 eta1 eta2 eta3 eta4 real
% new_state = [csi1; csi2; eta1; eta2; eta3; eta4];
% 
% syms v real
% u = -(Lf2H/LgLfH) + (1/LgLfH)*v;
% % u = -(Lf3H/Lg2LfH) + (1/Lg2LfH)*v;
% 
% % matlabFunction(u, 'File', 'alpha3rd', 'Vars', [theta; theta_dot; m; g; k; L; D; beta_r; beta; v], 'Outputs', {'alphaFL'}, 'Optimize', false)
% 
% %% Zero Dynamics
% disp("Zero Dynamics")
% f_zero = subs([eta3; eta4; subs([zeros(2, 4) eye(2, 2)]*F, x, inv(diffPhi)*new_state)], [csi1; csi2], zeros(2, 1));
% % Subs Parameters
% f_zero = subs(f_zero, [m; g; k; L; D; beta; beta_r], [1; 9.81; 1; 1; 0.1; 0.1; 0.5]);
% % eta3 = eta4 == 0
% f_zero = subs(f_zero, [eta3; eta4], [0; 0]);
% f_zero = f_zero(3:end);
% 
% equilibria_equation = simplify(f_zero) == zeros(2, 1);
% disp("Equilibria");
% n_try = 10;
% for j = 1:n_try
%     disp(num2str(j) + "-th try")
%     solutions = vpasolve(equilibria_equation, [eta1; eta2], 'Random', true);
%     if(isempty(solutions.eta1))
%         equilibria(j, 1) = nan;
%         equilibria(j, 2) = nan;
%         equilibria(j, 3) = nan;
%         equilibria(j, 4) = nan;
%     else
%         equilibria(j, 1) = double(solutions.eta1);
%         equilibria(j, 2) = double(solutions.eta2);
%         equilibria(j, 3) = double(solutions.eta3);
%         equilibria(j, 4) = double(solutions.eta4);
%     end
% end
% 
% save("equilibriaZeroAlpha.mat", "equilibria");

%% Cartesian Output
h = simplify(cos(theta_r)*subs(x_s, s, 1) - sin(theta_r)*subs(y_s, s, 1));
LfH = jacobian(h, x)*F;
LgH = jacobian(h, x)*G;
Lf2H = jacobian(LfH, x)*F;
LgLfH = jacobian(LfH, x)*G;

% syms v real
% u = -(Lf2H/LgLfH) + (1/LgLfH)*v;
% disp("Save control in a matlab Function");
% matlabFunction(u, 'File', 'xsFL', 'Vars', [theta; theta_dot; m; g; k; L; D; beta_r; beta; v], 'Outputs', {'xsFL'})
