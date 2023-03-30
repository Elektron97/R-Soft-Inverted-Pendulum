%%%%%%%%% Validazione Matrici Dinamiche %%%%%%%%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");
addpath(genpath("Della Santina"));
addpath("origin_soft_pendulum");
%% Parameters
L = 1;
D = 0.1;

m = 1;
g = 9.81;
beta = 0.1;
k = 1;

threshold = 1e-3;

%% Computational Cost Analysis 1D
% tic
% inertiaMatrix(pi/4, pi/4, -pi/4, m, L, D);
% inertia_comp = toc; %[s]
% disp("Inertia Matrix Computational time:" + num2str(inertia_comp) + " seconds.")
% 
% tic
% gravityVector(pi/4, pi/4, -pi/4, m, g, L, D);
% gravity_comp = toc; %[s]
% disp("Gravity Vector Computational time:" + num2str(gravity_comp) + " seconds.")

%% Computational Cost Analysis 2D
abs_range = 10;
theta1 = -abs_range:0.1:abs_range;
theta0 = zeros(size(theta1));
thetar = zeros(size(theta1));
n_try = 100;

for i=1:length(theta1)
    %Not Singular Configurations
    if(abs(theta1(i)) <= threshold)
       theta1_sing = threshold;
       
       if(theta0(i) <= threshold)
            theta0_sing = threshold;
       else
            theta0_sing = theta0(i);  
       end   
       
    else
        theta0_sing = theta0(i);
        theta1_sing = theta1(i);
    end

    % Several Tries to have a relevant statistic
    for j=1:n_try
        tic
        inertiaMatrix(thetar(i), theta0_sing, theta1_sing, m, L, D);
        realizations(j) = toc;
    end
    var_comp = std(realizations);
    inertia_comp(i) = mean(realizations); %[s]
end

%Plot
figure
plot(theta1, 1000*inertia_comp, 'LineWidth', 2)
% hold on
% plot(theta1, 1000*(inertia_comp + 3*var_comp), 'Color', [.7, .7, .7])
% plot(theta1, 1000*(inertia_comp - 3*var_comp), 'Color', [.7, .7, .7])
% % errorbar(theta1, 1000*inertia_comp, 1000*var_comp, 'Color', [0.7, 0.7, 0.7])
% hold off
grid on
xlabel("\theta_1 [rad]")
ylabel("Avg. Comp. Time [ms]")
title("Average Computational Time of Inertia Matrix")

%% Computational Cost Analysis 3D
% % Define of Range of Values
% abs_range = 10;
% % n_try = 100;
% n_try = 1;
% [THETA0, THETA1] = meshgrid(-abs_range:0.1:abs_range, -abs_range:0.1:abs_range);
% THETAR = 0*THETA0;
% 
% for i = 1:size(THETA0, 1)
%     for j = 1:size(THETA0, 2)
%         
%         % Not Singular Configurations
%         if(abs(THETA1(i, j)) <= threshold)
%            theta1_sing = threshold;
%            
%            if(abs(THETA0(i, j)) <= threshold)
%                 theta0_sing = threshold;
%            else
%                 theta0_sing = THETA0(i, j);  
%            end   
%            
%         else
%             theta0_sing = THETA0(i, j);
%             theta1_sing = THETA1(i, j);
%         end
%         % Several Tries to have a relevant statistic
%         for k=1:n_try
%             tic
%             inertiaMatrix(THETAR(i, j), theta0_sing, theta1_sing, m, L, D);
%             realizations(k) = toc;
%         end
%         inertia_comp(i, j) = mean(realizations); %[s]
%     end
% end
% 
% %% Plotting
% figure
% s = surf(THETA0, THETA1, inertia_comp);
% s.EdgeColor = 'none';
% xlabel("\theta_0 [rad]");
% ylabel("\theta_1 [rad]");
% zlabel("Avg. Comp. Time [s]");
% grid on
% title("Average Computational Time of Inertia Matrix")
% % xlim([-8, 8])
% % ylim([-8, 8])
% % zlim([0 0.01])

