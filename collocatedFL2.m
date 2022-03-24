%%%%%%%%%% Collocated Feedback %%%%%%%
% Output: alpha1
function tau_r = collocatedFL2(state, state_dot, Kd, m, g, k, L, D, beta_r, beta, Kp, alpha_des)
    
    thresh = 1e-3;
    %No Identically Zero
    for i = 1:length(state)
       if(abs(state(i)) < thresh)
           state(i) = thresh;
       end
    end

    %% Compute Dynamic Matrices
    B = inertiaMatrix(state(1), state(2), state(3), m, L, D);
    C = coriolisMatrix(state(1), state(2), state(3), state_dot(1), state_dot(2), state_dot(3), m, L, D);
    G = gravityVector(state(1), state(2), state(3), m, g, L, D);
    K = elasticMatrix(k);
    Damp = dampingMatrix(beta, beta_r);

    h = C*state_dot + K*state + Damp*state_dot;

    %% Block Matrices
    %Inertia
    Brr = B(1, 1);
    Bro = B(1, 2:3);
    Boo = B(2:3, 2:3);

    % Gravity
    Gr = G(1);
    Go = G(2:3);

    % Other terms
    hr = h(1);
    ho = h(2:3);

    %% Change of Variables
    inv_Boo = pinv(Boo);

    Brr_tilde = Brr - Bro*inv_Boo*(Bro');

    %% Block Inverse of B
    invBrr = pinv(Brr_tilde);
    invBro = -invBrr*(Bro*inv_Boo);
    invBor = -inv_Boo*(Bro')*(invBrr);
    invBoo = inv_Boo + inv_Boo*(Bro')*invBrr*Bro*inv_Boo;

    %% PD law
    csi1 = [1 1 1/2]*state;
    csi2 = [1 1 1/2]*state_dot;
    v = -Kd*csi2 + Kp*(alpha_des - csi1);

    %% Final Torque
    % Numerical Issues
%     tau_r = pinv(invBrr + [1 1/2]*invBor)*(v + [1 1 1/2]*(B\(h+G)));
    tau_r = pinv(invBrr + [1 1/2]*invBor)*(v + [1 1 1/2]*(pinv(B)*(h+G)));
end