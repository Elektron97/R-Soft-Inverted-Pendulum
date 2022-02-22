function ddx = Soft_dynamics(theta, theta_dot, tau, m, g, L, D, k, beta)

thresh = 1e-5;
%No Identically Zero
for i = 1:length(theta)
   if(abs(theta(i)) < thresh)
       theta(i) = thresh;
   end
end

%Origin Soft Inverted Pendulum
myB = originInertiaMatrix(theta(1), theta(2), m, L, D);
myG = originGravityVector(theta(1), theta(2), m, g, L, D);
myC = originCoriolisMatrix(theta(1), theta(2), theta_dot(1), theta_dot(2), m, L, D);
myD = originDampingMatrix(beta);
myK = originElasticMatrix(k);

ddx = -myB\(myC*theta_dot + myG + myK*theta + myD*theta_dot - henkelMatrix(2)*[1; 0]*tau);
end