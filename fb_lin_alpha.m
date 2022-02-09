function tau = fb_lin_alpha(t,dt)
%
% function tau = fb_lin_alpha(t,dt)
%
% This function implements a collocated fb linearization of the variable
% alpha(1) = theta_0 + theta_1/2
%
global beta k

if abs(t(1)) < 1e-3
    t(1) = 1e-3;
end

if abs(t(2)) < 1e-3
    t(2) = 1e-3;
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

%   b_bar = B_tilde(2,2) - (B_tilde(2,1)/B_tilde(1,1))*B_tilde(2,1);
%   h_bar =       h(2)   - (B_tilde(2,1)/B_tilde(1,1))*h(1);
% phi_bar =     phi(2)   - (B_tilde(2,1)/B_tilde(1,1))*phi(1);
  b_bar = B_tilde(1,1) - (B_tilde(2,1)/B_tilde(2,2))*B_tilde(2,1);
  h_bar =       h(1)   - (B_tilde(2,1)/B_tilde(2,2))*h(2);
phi_bar =     phi(1)   - (B_tilde(2,1)/B_tilde(2,2))*phi(2);

% g_P = 1;
% g_D = 1;
g_P = 10000;
g_D = 1000;
u = -g_P*q(1) -g_D*dq(1);

tau = b_bar*u + h_bar + phi_bar;

