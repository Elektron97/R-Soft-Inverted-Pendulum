%%%%%%%%%%%% Feedback Linearization %%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Model
% Rsoft_model

%% Redefine state variables
x1 = theta_r;
x2 = [theta0; theta1];
x3 = theta_r_dot;
x4 = [theta0_dot; theta1_dot];
x = [x1; x2; x3; x4];
%% Output
h = [1 1 1/2 0 0 0]*x;
