function tau_r = collocatedFL(state, state_dot, Kd, m, g, k, L, D, beta_r, beta, Kp, thetaR_des, thetaR_dot_des, thetaR_2dot_des)
    
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
    hr_tilde = hr - Bro*inv_Boo*ho;
    Gr_tilde = Gr - Bro*inv_Boo*Go;

    %% PD law
    v = thetaR_2dot_des + Kd*(thetaR_dot_des - state_dot(1)) + Kp*(thetaR_des - state(1));

    %% Final Torque
    tau_r = Brr_tilde*v + hr_tilde + Gr_tilde;
end