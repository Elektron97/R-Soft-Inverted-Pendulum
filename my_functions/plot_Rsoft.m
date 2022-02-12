%%%%%%%%%%%%%%%% Plot R-Soft Inverted Pendulum %%%%%%%%%%%%%%%
function plot_Rsoft(theta, L, D, length_arrow)
%% Local Variables
theta_r = theta(1);
theta0 = theta(2);
theta1 = theta(3);

%% Compute Forward Kinematics
Ri0 = my_rot(theta_r, 'z');
Ti0 = blkdiag(Ri0, 1);

%% Plot R-Soft Inverted Pendulum
%Inertial Frame
Ti = eye(4);
trplot(Ti, 'frame', 'I', 'color', 'k', 'length', length_arrow)
hold on

%S0 Frame
trplot(Ti0, 'frame', 'S0', 'length', length_arrow, 'color', 'r')

%Ss Frames
s_step = 0.01;
d_step = 0.01;

s = 0:s_step:1;

for i = 1:length(s)
    %% Affine Curvature
%     K = theta0 + theta1*s;

    %% Orientation of SoR {Ss} w.r.t. {S0}
    alpha(i) = theta0*s(i) + 0.5*theta1*(s(i)^2);

    if(theta1 == 0)
        if(theta0 == 0)
            %Limit results
            x_s = 0;
            y_s = L*s(i);
        else
            %Limit results
            x_s = L*(cos(s(i)*theta0)-1)/theta0;
            y_s = L*sin(s(i)*theta0)/theta0;
        end
       
    else

        fresn_sin1 = fresnels_approx( (theta0 + s(i)*theta1) *  sqrt(1/(pi*theta1)) );
        fresn_sin2 = fresnels_approx(theta0 *  sqrt(1/(pi*theta1)) );

        fresn_cos1 = fresnelc_approx((theta0 + s(i)*theta1) * sqrt(1/(pi*theta1)));
        fresn_cos2 = fresnelc_approx(theta0 * sqrt(1/(pi*theta1)) );

        x_s = L*( sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) - cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));

        y_s = L*(cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) + sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));

    end
    p_s_hom = [x_s; y_s; 0; 1];
    p_sI_hom = Ti0*p_s_hom;
    p_s(:, i) = p_sI_hom(1:3);

    T0s = [my_rot(alpha(i), 'z') p_s_hom(1:3); zeros(1, 3) 1] ;
    
    %% Plot Ss frame
    if(s(i)==1)
        trplot(Ti0*T0s, 'frame', 'S1', 'length', length_arrow, 'color', 'r')
    end
    
    for d=-0.5:d_step:0.5
        p_sd_hom = Ti0*T0s*[d*D; 0; 0; 1];
        
        if(d == -0.5)
            p_sd_down(:, i) = p_sd_hom(1:3);
        elseif(d == 0.5)
            p_sd_up(:, i) = p_sd_hom(1:3);
        end
    end
end

% 1D Soft Segment
plot3(p_s(1, :), p_s(2, :), p_s(3, :), 'linewidth', 2, 'color', [0, 0.4470, 0.7410])

% Thickness
plot3(p_sd_down(1, :), p_sd_down(2, :), p_sd_down(3, :), 'color', [0, 0.4470, 0.7410])
plot3(p_sd_up(1, :), p_sd_up(2, :), p_sd_up(3, :), 'color', [0, 0.4470, 0.7410])
hold off

xlim([-1 1])
ylim([-1 1])
grid on
end