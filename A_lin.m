function A = A_lin(theta, m, g, k, L, D, beta, beta_r)

thresh = 1e-3;
%No Identically Zero
for i = 1:length(theta)
   if(abs(theta(i)) < thresh)
       theta(i) = thresh;
   end
end

%R-Soft Inverted Pendulum Linearized
B = inertiaMatrix(theta(1), theta(2), theta(3), m, L, D);
Damp = dampingMatrix(beta, beta_r);
Stiff_Mat = stiffMatrix(theta(1), theta(2), theta(3), m, g, k, L, D);


A = [zeros(3, 3), eye(3,3); ...
    -B\Stiff_Mat, -B\Damp];
end