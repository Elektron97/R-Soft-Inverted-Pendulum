%%%%%%%%%% Equilibria Analysis %%%%%%%%%%
clear all
close all
clc

%% AddPath functions
addpath("my_functions");

%% Load Equilibria
load("equilibria.mat");

%% Filter Equilibria
for i = 1:length(equilibria)
   filtered_equilibria{i} = filterEquilibria(equilibria{i});
end

%% Stability Properties
tau_r = -5:0.1:5;

for i = 1:length(tau_r)
    for j = 1:size(filtered_equilibria{i}, 1)     
        stability{i}(j) = isHurwitz(A_lin(filtered_equilibria{i}(j, :)', 1, 9.81, 1, 1, 0.1, 0.1, 0.5)); 
%         stability{i}(j) = isHurwitz(A_lin2(filtered_equilibria{i}(j, :)', tau_r(i), 1, 9.81, 1, 1, 0.1, 0.1, 0.5)); 
    end
    disp("Stabilita' per tau_r: " + num2str(tau_r(i)))
end
%% Plot Equilibria
f1 = figure;
f2 = figure;
f3 = figure;

for i = 1:length(tau_r)
    figure(f1)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)  
        
        switch(stability{i}(j))
            case -1
                plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'rx')
                
            case 0
                plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'gx')
                
            case 1
                plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'bx')
        end 
    end
%     grid on
%     title("Equilibria: \theta_r")
%     xlabel("\tau")
%     ylabel("\theta_r [rad]")
%     hold off
%     legend("Instable", "Stable")

    figure(f2)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        switch(stability{i}(j))
            case -1
                plot(tau_r(i), filtered_equilibria{i}(j, 2), 'rx')
                
            case 0
                plot(tau_r(i), filtered_equilibria{i}(j, 2), 'gx')
                
            case 1
                plot(tau_r(i), filtered_equilibria{i}(j, 2), 'bx')
        end 
        
    end
%     grid on
%     title("Equilibria: \theta_0")
%     xlabel("\tau")
%     ylabel("\theta_0 [rad]")
%     hold off
%     legend("Instable", "Stable")

    figure(f3)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        switch(stability{i}(j))
            case -1
                plot(tau_r(i), filtered_equilibria{i}(j, 3), 'rx')
                
            case 0
                plot(tau_r(i), filtered_equilibria{i}(j, 3), 'gx')
                
            case 1
                plot(tau_r(i), filtered_equilibria{i}(j, 3), 'bx')
        end 
        
    end
%     grid on
%     title("Equilibria: \theta_1")
%     xlabel("\tau")
%     ylabel("\theta_1 [rad]")
%     hold off
%     legend("Unstable", "Stable")
end

figure(f1)
grid on
title("Equilibria: \theta_r")
xlabel("\tau")
ylabel("\theta_r [rad]")
hold off
legend("Unstable", "Stable")

figure(f2)
grid on
title("Equilibria: \theta_0")
xlabel("\tau")
ylabel("\theta_0 [rad]")
hold off
legend("Unstable", "Stable")

figure(f3)
grid on
title("Equilibria: \theta_1")
xlabel("\tau")
ylabel("\theta_1 [rad]")
hold off
legend("Unstable", "Stable")

%% Controllability
for i = 1:length(tau_r)
    for j = 1:size(filtered_equilibria{i}, 1)     
        A = A_lin(filtered_equilibria{i}(j, :)', 1, 9.81, 1, 1, 0.1, 0.1, 0.5);
        B = B_lin(filtered_equilibria{i}(j, :)', [1; 0; 0], 1, 1, 0.1);
        
        local_control{i}(j) = det(ctrb(A, B)); 
    end
    local_control{i}(local_control{i} ~= 0) = true;
    local_control{i}(local_control{i} == 0) = false;
end

%% Plot Only Stable and Controllable Equilibria
f4 = figure;
f5 = figure;
f6 = figure;

for i = 1:length(tau_r)
    figure(f4)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)  
        if((stability{i}(j) == 1) && local_control{i}(j))
            plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'gx')
        end
    end

    figure(f5)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        if((stability{i}(j) == 1) && local_control{i}(j))
            plot(tau_r(i), filtered_equilibria{i}(j, 2), 'gx')
        end
    end

    figure(f6)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        if((stability{i}(j) == 1) && local_control{i}(j))
            plot(tau_r(i), filtered_equilibria{i}(j, 3), 'gx')
        end
    end
end

figure(f4)
grid on
title("Equilibria: \theta_r")
xlabel("\tau")
ylabel("\theta_r [rad]")
hold off
legend("Stable and Controllable")

figure(f5)
grid on
title("Equilibria: \theta_0")
xlabel("\tau")
ylabel("\theta_0 [rad]")
hold off
legend("Stable and Controllable")

figure(f6)
grid on
title("Equilibria: \theta_1")
xlabel("\tau")
ylabel("\theta_1 [rad]")
hold off
legend("Stable and Controllable")

%% Observability
for i = 1:length(tau_r)
    for j = 1:size(filtered_equilibria{i}, 1)     
        A = A_lin(filtered_equilibria{i}(j, :)', 1, 9.81, 1, 1, 0.1, 0.1, 0.5);
        C = C_lin(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), filtered_equilibria{i}(j, 3), 1);
        C = [[1 0 0]; C];
        rankObs = rank(obsv(A, [C, zeros(3, 3)]));

        if(rankObs == 6)
            local_obs{i}(j) = true;  
        else
            local_obs{i}(j) = false;
        end
    end
end

%% Plot Observability
f7 = figure;
f8 = figure;
f9 = figure;

for i = 1:length(tau_r)
    figure(f7)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)  
        if(local_obs{i}(j))
            plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'mx')
        else
            plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'cx')
        end
    end

    figure(f8)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        if(local_obs{i}(j))
            plot(tau_r(i), filtered_equilibria{i}(j, 2), 'mx')
        else
            plot(tau_r(i), filtered_equilibria{i}(j, 2), 'cx')
        end
    end

    figure(f9)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        if(local_obs{i}(j))
            plot(tau_r(i), filtered_equilibria{i}(j, 3), 'mx')
        else
            plot(tau_r(i), filtered_equilibria{i}(j, 3), 'cx')
        end
    end
end

figure(f7)
grid on
title("Equilibria: \theta_r")
xlabel("\tau")
ylabel("\theta_r [rad]")
hold off
legend("Observable")

figure(f8)
grid on
title("Equilibria: \theta_0")
xlabel("\tau")
ylabel("\theta_0 [rad]")
hold off
legend("Observable")

figure(f9)
grid on
title("Equilibria: \theta_1")
xlabel("\tau")
ylabel("\theta_1 [rad]")
hold off
legend("Observable")

%% Rigid
% load("equilibria_elastic.mat");
% 
% %% Filter Equilibria
% theta_r_eq = [];
% theta0_eq = [];
% theta1_eq = [];
% 
% for i = 1:length(equilibria)
%    filtered_equilibria{i} = filterEquilibria(equilibria{i});
% end
% 
% rigid = 1e-5:0.5:10;
% %% Plot Equilibria
% f1 = figure;
% f2 = figure;
% f3 = figure;
% 
% for i = 1:length(rigid)
%     figure(f1)
%     hold on
%     for j = 1:size(filtered_equilibria{i}, 1)
%         plot(rigid(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'rx')
%     end
%     grid on
%     title("Equilibria: \theta_r")
%     xlabel("k")
%     ylabel("\theta_r [rad]")
%     hold off
% 
%     figure(f2)
%     hold on
%     for j = 1:size(filtered_equilibria{i}, 1)
%         plot(rigid(i), wrapToPi(filtered_equilibria{i}(j, 2)), 'rx')
%     end
%     grid on
%     title("Equilibria: \theta_0")
%     xlabel("k")
%     ylabel("\theta_0 [rad]")
%     hold off
% 
%     figure(f3)
%     hold on
%     for j = 1:size(filtered_equilibria{i}, 1)
%         plot(rigid(i), wrapToPi(filtered_equilibria{i}(j, 3)), 'rx')
%     end
%     grid on
%     title("Equilibria: \theta_1")
%     xlabel("k")
%     ylabel("\theta_1 [rad]")
%     hold off
% end

%% Filter Equilibria
function equil = filterEquilibria(equilibria)
    equil = zeros(1, 3);
    
    for i = 1:size(equilibria, 1)
        if(equilibria(i, :) ~= equil)
            equil = [equil; equilibria(i, :)];
        end
    end
    
    equil = equil(2:end, :);
end