function ddx = RSoft_dynamics(theta, theta_dot, tau_r, m, g, L, D, k, beta)

thresh = 1e-5;
%No Identically Zero
for i = 1:length(theta)
   if(abs(theta(i)) < thresh)
       theta(i) = thresh;
   end
end

%Origin Soft Inverted Pendulum
B = inertiaMatrix(theta(1), theta(2), theta(3), m, L, D);
G = gravityVector(theta(1), theta(2), theta(3), m, g, L, D);
C = coriolisMatrix(theta(1), theta(2), theta(3), theta_dot(1), theta_dot(2), theta_dot(3), m, L, D);
Damp = dampingMatrix(beta);
K = elasticMatrix(k);

ddx = -B\(C*theta_dot + G + K*theta + Damp*theta_dot - [1; 0; 0]*tau_r);
end