%%%%%%%%%%%%%%%% Pick and Place planner %%%%%%%%%%%%%
clear all
close all
clc

%% Addpath functions
addpath("my_functions");

%% Load Equilibria of Controlled System
load("equilibria_phi1.mat");

%% Filter Equilibria
phi = 0:0.1:2*pi;
for i = 1:length(equilibria)
   filtered_equilibria{i} = filterEquilibria(equilibria{i});
end

%% Stability
for i = 1:length(phi)
    for j = 1:size(filtered_equilibria{i}, 1)     
        stability{i}(j) = isPositiveDef(originStiffMat(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), 1, 9.81, 1, 1, 0.1, phi(i))); 
    end
end

%% Forward Kinematics
for i = 1:length(phi)
    for j = 1:size(filtered_equilibria{i}, 1)
        if stability{i}(j)
            workSpace{i}(j, :) = fwdKinRSIP([phi(i); filtered_equilibria{i}(j, :)'], 1, 0, 1, 0.1);
        else
            workSpace{i}(j, :) = nan*ones(3, 1);
        end
    end
end

%% Work Space
figure
hold on
for i = 1:length(phi)
    for j = 1:size(filtered_equilibria{i}, 1)
        if stability{i}(j)
            plot(workSpace{i}(j, 1), workSpace{i}(j, 2), 'gx')
        end
    end
end
hold off
grid on
xlabel("x [m]")
ylabel("y [m]")
title("Workspace fo R-sip")




%% Functions
function equil = filterEquilibria(equilibria)
    equil = zeros(1, 2);
    
    for i = 1:size(equilibria, 1)
        if(equilibria(i, :) ~= equil)
            equil = [equil; equilibria(i, :)];
        end
    end
    
    equil = equil(2:end, :);
end