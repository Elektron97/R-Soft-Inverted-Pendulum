%%% Object Analysis %%%
clear all
close all
clc

%% AddPath functions
addpath("my_functions");

%% Load Equilibria
% Without Object
load("equilibria_phi1.mat");
% Reassign
equilibria0 = equilibria;

% 0.1 [Kg]
load("equilibria_phi_with_obj.mat");
% Reassign
equilibria1 = equilibria;

% 0.3 [Kg]
load("equilibria_phi_with_obj_3.mat");
% Reassign
equilibria3 = equilibria;

clear equilibria

equilibria_obj = {equilibria0, equilibria1, equilibria3};

%% Filter Equilibria
% Range of phi
phi = 0:0.1:2*pi;
m_objs = [0.0, 0.1, 0.3];

% Cycle on the different mass values
for k = 1:length(equilibria_obj)
    for i = 1:length(equilibria_obj{k})
       filtered_equilibria{k}{i} = filterEquilibria(equilibria_obj{k}{i});
    end


%% Stability
    for i = 1:length(phi)
        for j = 1:size(filtered_equilibria{k}{i}, 1)     
            stability{k}{i}(j) = isPositiveDef(originStiffMat(filtered_equilibria{k}{i}(j, 1), filtered_equilibria{k}{i}(j, 2), 1 + m_objs(k), 9.81, 1, 1, 0.1, phi(i))); 
        end
    end
end

%% Plot WrapToPi equilibria
% Color array
stable_colors = ["#0072BD", "#00A4D9", "#00D1D8", "#75FAC8"];
unstable_colors = ["#E3170A", "#FF5358", "#FF889B", "#FFBCDA"];

f3 = figure;
f4 = figure;

for k = 1:length(equilibria_obj)
    for i = 1:length(phi)/2
        figure(f3)
        hold on
        for j = 1:size(filtered_equilibria{k}{i}, 1)
            
            if(stability{k}{i}(j))
                plot(phi(i), filtered_equilibria{k}{i}(j, 1), 'o', 'Color', stable_colors(k), 'LineWidth', 2.0)
            else
                plot(phi(i), filtered_equilibria{k}{i}(j, 1), 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0)
            end 
            
        end
    
        figure(f4)
        hold on
        for j = 1:size(filtered_equilibria{k}{i}, 1)
            
            if(stability{k}{i}(j))
                plot(phi(i), filtered_equilibria{k}{i}(j, 2), 'o', 'Color', stable_colors(k), 'LineWidth', 2.0)
            else
                plot(phi(i), filtered_equilibria{k}{i}(j, 2), 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0)
            end 
            
        end
    end
    
    for i = 32:length(phi)
        figure(f3)
        hold on
        for j = 1:size(filtered_equilibria{k}{i}, 1)
            
            if(stability{k}{i}(j))
                plot(phi(i) - 2*pi, filtered_equilibria{k}{i}(j, 1), 'o', 'Color', stable_colors(k), 'LineWidth', 2.0)
            else
                plot(phi(i) - 2*pi, filtered_equilibria{k}{i}(j, 1), 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0)
            end 
            
        end
    
        figure(f4)
        hold on
        for j = 1:size(filtered_equilibria{k}{i}, 1)
            
            if(stability{k}{i}(j))
                plot(phi(i) - 2*pi, filtered_equilibria{k}{i}(j, 2), 'o', 'Color', stable_colors(k), 'LineWidth', 2.0)
            else
                plot(phi(i) - 2*pi, filtered_equilibria{k}{i}(j, 2), 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0)
            end 
            
        end
    end

    % Dummy plots to create legends
    figure(f3)
    hold on
        plt_stable0{k} = plot(nan, 'o', 'Color', stable_colors(k), 'LineWidth', 2.0);
        plt_unstable0{k} = plot(nan, 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0);
    hold off

    figure(f4)
    hold on
        plt_stable1{k} = plot(nan, 'o', 'Color', stable_colors(k), 'LineWidth', 2.0);
        plt_unstable1{k} = plot(nan, 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0);
    hold off
end

figure(f3)
grid on
legend([plt_stable0{:}, plt_unstable0{:}], ...
        {'Stable 0.0 [Kg]', 'Stable 0.1 [Kg]', 'Stable 0.3 [Kg]', ...
         'Unstable 0.0 [Kg]', 'Unstable 0.1 [Kg]', 'Unstable 0.3 [Kg]'}, 'NumColumns', 2);
xlabel("\theta_{r, d} [rad]")
ylabel("\theta_0 [rad]")
xlim([-pi, pi])


figure(f4)
grid on
xlabel("\theta_{r, d} [rad]")
ylabel("\theta_1 [rad]")
xlim([-pi, pi])
legend([plt_stable1{:}, plt_unstable1{:}], ...
                    {'Stable 0.0 [Kg]', 'Stable 0.1 [Kg]', 'Stable 0.3 [Kg]', ...
                     'Unstable 0.0 [Kg]', 'Unstable 0.1 [Kg]', 'Unstable 0.3 [Kg]'}, 'NumColumns', 2);

%% 3D Equilibria with phase space trajectory
f5 = figure;

for k = [1 3]
    for i = 1:length(phi)/2
        figure(f5)
        hold on
        for j = 1:size(filtered_equilibria{k}{i}, 1)
    
            if(stability{k}{i}(j))
                plot3(filtered_equilibria{k}{i}(j, 1), filtered_equilibria{k}{i}(j, 2), phi(i), 'o', 'Color', stable_colors(k), 'LineWidth', 1.5)
            else
                plot3(filtered_equilibria{k}{i}(j, 1), filtered_equilibria{k}{i}(j, 2), phi(i), 'x', 'Color', unstable_colors(k), 'LineWidth', 1.5)
            end 
    
        end
    end
    
    for i = 32:length(phi)
        figure(f5)
        hold on
        for j = 1:size(filtered_equilibria{k}{i}, 1)
    
            if(stability{k}{i}(j))
                plot3(filtered_equilibria{k}{i}(j, 1), filtered_equilibria{k}{i}(j, 2), phi(i)-2*pi, 'o', 'Color', stable_colors(k), 'LineWidth', 1.5)
            else
                plot3(filtered_equilibria{k}{i}(j, 1), filtered_equilibria{k}{i}(j, 2), phi(i)-2*pi, 'x', 'Color', unstable_colors(k), 'LineWidth', 1.5)
            end 
    
        end
    end

    hold on
        plt_stable2{k} = plot(nan, 'o', 'Color', stable_colors(k), 'LineWidth', 2.0);
        plt_unstable2{k} = plot(nan, 'x', 'Color', unstable_colors(k), 'LineWidth', 2.0);
end
load("result_PP_with_obj.mat");
plot3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1), ...
      'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5)
% quiver3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1), ...
%         result.simout3.data(:, 2), result.simout3.data(:, 3), result.simout3.data(:, 1),'Color', [0.4940 0.1840 0.5560])

% % Start and Final Points
% Start
marker_width = 2.5;
marker_size = 11;

plot3(result.simout.data(1, 2), result.simout.data(1, 3), result.simout.data(1, 1), ...
        'o', 'LineWidth', marker_width, 'MarkerSize', marker_size, 'Color', "#77AC30")
% End
plot3(result.simout.data(end, 2), result.simout.data(end, 3), result.simout.data(end, 1), ...
        '^', 'LineWidth', marker_width, 'MarkerSize', marker_size, 'Color', "#7E2F8E")
hold off

figure(f5)
grid on
xlabel("\theta_0 [rad]")
ylabel("\theta_1 [rad]")
zlabel("\theta_{rd} [rad]")
view(109, 15)
legend([plt_stable2{[1, 3]}, plt_unstable2{[1, 3]}], ...
        {'Stable 0.0 [Kg]', 'Stable 0.3 [Kg]', ...
         'Unstable 0.0 [Kg]', 'Unstable 0.3 [Kg]'}, 'NumColumns', 2);
zlim([-pi pi])

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