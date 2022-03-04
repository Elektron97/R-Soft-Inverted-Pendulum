%%%%% Test Linearization: Derivation of Inverse of Inertia Matrix %%%%%
inv_B = simplify(inv(B));
disp("Inversa di B calcolata");
deriv_inv_B = diffMatrix(inv_B, theta);
disp("Derivata di B calcolata");

