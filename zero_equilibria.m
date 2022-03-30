%%%%%%%%%% Zero Dynamics Equilibria Analysis %%%%%%%%%%
clear all
close all
clc

%% AddPath functions
addpath("my_functions");

%% Load Equilibria
k = 4;

if k == 1
    load("equilibria_phi1.mat");    
elseif k == 4
    load("equilibria_phi4.mat");
end
%% Filter Equilibria
for i = 1:length(equilibria)
   filtered_equilibria{i} = filterEquilibria(equilibria{i});
end

%% Plot Equilibria
phi = 0:0.1:2*pi;
% f1 = figure;
% f2 = figure;
% 
% for i = 1:length(phi)
%     figure(f1)
%     hold on
%     plot(phi(i), filtered_equilibria{i}(:, 1), 'rx');
%     hold off
%     
%     figure(f2)
%     hold on
%     plot(phi(i), equilibria{i}(:, 2), 'rx');
%     hold off
% end
% 
% figure(f1)
% grid on
% xlabel("\theta_{r, d} [rad]")
% ylabel("\theta_0 [rad]")
% title("Equilibria of Zero Dynamics: \theta_0")
% 
% figure(f2)
% grid on
% xlabel("\theta_{r, d} [rad]")
% ylabel("\theta_1 [rad]")
% title("Equilibria of Zero Dynamics: \theta_1")

%% Stability
for i = 1:length(phi)
    for j = 1:size(filtered_equilibria{i}, 1)     
        stability{i}(j) = isPositiveDef(originStiffMat(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), 1, 9.81, k, 1, 0.1, phi(i))); 
    end
end

%% Plot Equilibria
f1 = figure;
f2 = figure;

for i = 1:length(phi)
    figure(f1)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot(phi(i), filtered_equilibria{i}(j, 1), 'bx')
        else
            plot(phi(i), filtered_equilibria{i}(j, 1), 'rx')
        end 
        
    end

    figure(f2)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot(phi(i), filtered_equilibria{i}(j, 2), 'bx')
        else
            plot(phi(i), filtered_equilibria{i}(j, 2), 'rx')
        end 
        
    end
end

figure(f1)
grid on
xlabel("\theta_{r, d} [rad]")
ylabel("\theta_0 [rad]")
% legend("Unstable", "Stable")
xlim([phi(1), phi(end)])
title("Equilibria of Zero Dynamics: \theta_0")

figure(f2)
grid on
xlabel("\theta_{r, d} [rad]")
ylabel("\theta_1 [rad]")
% legend("Unstable", "Stable")
xlim([phi(1), phi(end)])
title("Equilibria of Zero Dynamics: \theta_1")

%% Function
function equil = filterEquilibria(equilibria)
    equil = zeros(1, 2);
    
    for i = 1:size(equilibria, 1)
        if(equilibria(i, :) ~= equil)
            equil = [equil; equilibria(i, :)];
        end
    end
    
    equil = equil(2:end, :);
end