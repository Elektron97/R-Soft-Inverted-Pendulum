%%%%%%%%%%%%%%%%% Regressor %%%%%%%%%%%%%%%%%
function Y = regressorSoftInverted(theta0, theta0_dot, theta0_des, theta0_dot_des, theta0_2dot_des, Lambda)
%% q0,r
theta0_r_dot = theta0_dot_des + Lambda*(theta0_des - theta0);
theta0_r_2dot = theta0_2dot_des + Lambda*(theta0_dot_des - theta0_dot);

%% Useful Variables
threshold = 1e-5;
if(abs(theta0) < threshold)
    sinc0 = 1;
    sinc1 = 0;
    sinc2 = -1/3;
else
    sinc0 = sin(theta0/2)/(theta0/2);
    sinc1 = cos(theta0/2)/theta0 - (2*sin(theta0/2))/(theta0^2);
    sinc2 = (4*sin(theta0/2))/(theta0^3) - sin(theta0/2)/(2*theta0) - (2*cos(theta0/2))/(theta0^2);
end
%% Regressor
y1 = (theta0_r_2dot/16)*(sinc1^2 + 4*sinc0^2 + 16) - (theta0_r_dot*theta0_dot/11)*sinc1*(sinc2 + 2*sinc0);
y2 = (-1/4)*(2*sinc0*cos(theta0/2) + sinc1*sin(theta0/2));
y3 = theta0;
y4 = theta0_r_dot;

Y = [y1, y2, y3, y4];
end