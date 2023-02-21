%%%%%%%%%% Equilibria Analysis for 2 Stiffness Values %%%%%%%%
clear all
close all
clc

%% AddPath functions
addpath("my_functions");

%% Load Equilibria
load("equilibria_stiffness.mat");

%% Filtering Equilibria
for i = 1:size(equilibria, 1)
    for j = 1:size(equilibria, 2)
        filtered_equilibria{i, j} = filterEquilibria(equilibria{i});
    end
end

%% Plot Equilibria varying Stiffness
% Stiffness Parameters
step_k = 0.5;
k_values = 1e-5:step_k:10;

f1 = figure;
f2 = figure;
f3 = figure;

for i = 1:size(equilibria, 1)
    for j = 1:size(equilibria, 2)
        % \theta_r
        figure(f1)
        hold on
        for h = 1:size(filtered_equilibria{i, j}, 1)  
             plot3(k_values(i), k_values(j), wrapToPi(filtered_equilibria{i, j}(h, 1)), 'rx')
        end

        % \theta_0
        figure(f2)
        hold on
        for h = 1:size(filtered_equilibria{i, j}, 1)  
             plot3(k_values(i), k_values(j), filtered_equilibria{i, j}(h, 2), 'rx')
        end

        % \theta_1
        figure(f3)
        hold on
        for h = 1:size(filtered_equilibria{i, j}, 1)  
             plot3(k_values(i), k_values(j), filtered_equilibria{i, j}(h, 3), 'rx')
        end
    end
end

figure(f1)
grid on
xlabel("Joint Stiffness k_r");
ylabel("Bending Stiffness k");
zlabel("Equilibria \theta_r");

figure(f2)
grid on
xlabel("Joint Stiffness k_r");
ylabel("Bending Stiffness k");
zlabel("Equilibria \theta_0");

figure(f3)
grid on
xlabel("Joint Stiffness k_r");
ylabel("Bending Stiffness k");
zlabel("Equilibria \theta_1");
%% Filter Equilibria Function
function equil = filterEquilibria(equilibria)
    equil = zeros(1, 3);
    
    for i = 1:size(equilibria, 1)
        if(equilibria(i, :) ~= equil)
            equil = [equil; equilibria(i, :)];
        end
    end
    
    equil = equil(2:end, :);
end