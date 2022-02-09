function tau = fb_lin_I_alpha(t,dt)
%
% function tau = fb_lin_I_alpha(t,dt)
%
% This function implements a collocated fb linearization of the variable
% \int_0^1 alpha(1) = theta_0/2 + theta_1/3
%
global beta k p_vect

if abs(t(1)) < 1e-5
    t(1) = sign(t(1))*1e-5;
end

if abs(t(2)) < 1e-5
    t(2) = sign(t(2))*1e-5;
end

 H = [1, 1/2; 1/2, 1/3];
iH = inv(H);
 q = H* t;
dq = H*dt;

B_tilde = iH*B_fcn(t)*iH;
C_tilde = iH*C_fcn(t,dt)*dt;
G_tilde = iH*G_fcn(t);
D_tilde = beta*dt;
K_tilde = k*t;

  h = C_tilde + D_tilde;
phi = G_tilde + K_tilde;

% B_tilde = iH*B_fcn(iH*q)*iH;
% C_tilde = 0*iH*C_fcn(iH*q,iH*dq)*iH;
% G_tilde = 0*iH*G_fcn(iH*q);
% K_tilde = k*iH*q;
% D_tilde = beta*iH*dq;
% 
% h = C_tilde*dq + D_tilde;
% phi = G_tilde + K_tilde;

  b_bar = B_tilde(2,1) - (B_tilde(1,1)/B_tilde(1,2))*B_tilde(2,2);
  h_bar =       h(1)   - (B_tilde(1,1)/B_tilde(1,2))*h(2);
phi_bar =     phi(1)   - (B_tilde(1,1)/B_tilde(1,2))*phi(2);

% g_D_1 = 1e4;%
% g_D_2 = 1e1;%
% g_P_1 = 1e3;%
% g_P_2 = 1e0;%
% g_D_1 = -p_vect(1);%
% g_D_2 = -p_vect(2);%
% g_P_1 = -p_vect(3);%
% g_P_2 = -p_vect(4);%
% u = 0.1; %-g_P_1*q(1) -g_D_1*dq(1) -g_P_2*q(2) -g_D_2*dq(2);
%p_vect = [1e3,  1e0, 1e4, 1e1];
u = -p_vect*[q;dq];

tau = b_bar*u + h_bar + phi_bar;
