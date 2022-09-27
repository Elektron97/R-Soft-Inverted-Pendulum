function pose2D = fwdKinRSIP(theta, s, d, L, D)

%% Local Variables
theta_r = theta(1);
theta0 = theta(2);
theta1 = theta(3);

%% Compute Forward Kinematics
Ri0 = my_rot(theta_r, 'z');
Ti0 = blkdiag(Ri0, 1);

    %% Orientation of SoR {Ss} w.r.t. {S0}
    alpha = theta0*s + 0.5*theta1*(s^2);

    if(theta1 == 0)
        if(theta0 == 0)
            %Limit results
            x_s = 0;
            y_s = L*s;
        else
            %Limit results
            x_s = L*(cos(s*theta0)-1)/theta0;
            y_s = L*sin(s*theta0)/theta0;
        end
       
    else

        fresn_sin1 = fresnels_approx( (theta0 + s*theta1) *  sqrt(1/(pi*theta1)) );
        fresn_sin2 = fresnels_approx(theta0 *  sqrt(1/(pi*theta1)) );

        fresn_cos1 = fresnelc_approx((theta0 + s*theta1) * sqrt(1/(pi*theta1)));
        fresn_cos2 = fresnelc_approx(theta0 * sqrt(1/(pi*theta1)) );

        x_s = L*( sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) - cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));

        y_s = L*(cos((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_cos1 - fresn_cos2) + sin((theta0^2)/(2*theta1)) * sqrt(pi/theta1) * (fresn_sin1 - fresn_sin2));

    end
    p_s_hom = [x_s; y_s; 0; 1];
    p_sI_hom = Ti0*p_s_hom;
    p_s = p_sI_hom(1:3);

    T0s = [my_rot(alpha, 'z') p_s_hom(1:3); zeros(1, 3) 1] ;
    p_sd_hom = Ti0*T0s*[d*D; 0; 0; 1];

    x_sd = p_sd_hom(1);
    y_sd = p_sd_hom(2);
    alpha_s = alpha + theta_r;

    pose2D = [x_sd; y_sd; alpha_s];

end