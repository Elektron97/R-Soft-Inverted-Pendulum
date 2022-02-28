%%%%%%%%%% Equilibria Analysis %%%%%%%%%%
clear all
close all
clc

%% Load Equilibria
load("equilibria.mat");

%% Filter Equilibria
for i = 1:length(equilibria)
   filtered_equilibria{i} = filterEquilibria(equilibria{i});
end

%% Plot Equilibria
tau_r = -5:0.1:5;

f1 = figure;
f2 = figure;
f3 = figure;

for i = 1:length(tau_r)
    figure(f1)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        plot(tau_r(i), wrapToPi(filtered_equilibria{i}(j, 1)), 'rx')
    end
    grid on
    title("Equilibria: \theta_r")
    xlabel("\tau")
    ylabel("\theta_r [rad]")
    hold off

    figure(f2)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        plot(tau_r(i), filtered_equilibria{i}(j, 2), 'rx')
    end
    grid on
    title("Equilibria: \theta_0")
    xlabel("\tau")
    ylabel("\theta_0 [rad]")
    hold off

    figure(f3)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        plot(tau_r(i), filtered_equilibria{i}(j, 3), 'rx')
    end
    grid on
    title("Equilibria: \theta_1")
    xlabel("\tau")
    ylabel("\theta_1 [rad]")
    hold off
end

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
% %% Plot Equilibria
% rigid = 1e-5:0.5:10;
% 
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