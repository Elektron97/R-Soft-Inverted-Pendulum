function B = B_lin(theta, S, m, L, D)

thresh = 1e-3;
%No Identically Zero
for i = 1:length(theta)
   if(abs(theta(i)) < thresh)
       theta(i) = thresh;
   end
end

%R-Soft Inverted Pendulum Linearized
M = inertiaMatrix(theta(1), theta(2), theta(3), m, L, D);


B = [zeros(3, 1); ...
    M\S];
end