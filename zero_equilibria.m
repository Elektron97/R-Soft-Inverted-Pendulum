%%%%%%%%%% Zero Dynamics Equilibria Analysis %%%%%%%%%%
clear all
close all
clc

%% AddPath functions
addpath("my_functions");

%% Load Equilibria
k = 1;

if k == 1
    load("equilibria_phi1.mat");    
elseif k == 4
    load("equilibria_phi4.mat");
end
%% Filter Equilibria
for i = 1:length(equilibria)
   filtered_equilibria{i} = filterEquilibria(equilibria{i});
end

phi = 0:0.1:2*pi;
%% Stability
for i = 1:length(phi)
    for j = 1:size(filtered_equilibria{i}, 1)     
        stability{i}(j) = isPositiveDef(originStiffMat(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), 1, 9.81, k, 1, 0.1, phi(i))); 
    end
end

%% Plot Equilibria
% f1 = figure;
% f2 = figure;
% 
% for i = 1:length(phi)
%     figure(f1)
%     hold on
%     for j = 1:size(filtered_equilibria{i}, 1)
%         
%         if(stability{i}(j))
%             plot(phi(i), filtered_equilibria{i}(j, 1), 'bx')
%         else
%             plot(phi(i), filtered_equilibria{i}(j, 1), 'rx')
%         end 
%         
%     end
% 
%     figure(f2)
%     hold on
%     for j = 1:size(filtered_equilibria{i}, 1)
%         
%         if(stability{i}(j))
%             plot(phi(i), filtered_equilibria{i}(j, 2), 'bx')
%         else
%             plot(phi(i), filtered_equilibria{i}(j, 2), 'rx')
%         end 
%         
%     end
% end
% 
% figure(f1)
% grid on
% xlabel("\theta_{r, d} [rad]")
% ylabel("\theta_0 [rad]")
% % legend("Unstable", "Stable")
% xlim([phi(1), phi(end)])
% title("Equilibria of Zero Dynamics: \theta_0")
% 
% figure(f2)
% grid on
% xlabel("\theta_{r, d} [rad]")
% ylabel("\theta_1 [rad]")
% % legend("Unstable", "Stable")
% xlim([phi(1), phi(end)])
% title("Equilibria of Zero Dynamics: \theta_1")

%% WrapToPi equilibria
f3 = figure;
f4 = figure;

for i = 1:length(phi)/2
    figure(f3)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot(phi(i), filtered_equilibria{i}(j, 1), 'bx')
        else
            plot(phi(i), filtered_equilibria{i}(j, 1), 'rx')
        end 
        
    end

    figure(f4)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot(phi(i), filtered_equilibria{i}(j, 2), 'bx')
        else
            plot(phi(i), filtered_equilibria{i}(j, 2), 'rx')
        end 
        
    end
end

for i = 32:length(phi)
    figure(f3)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot(phi(i) - 2*pi, filtered_equilibria{i}(j, 1), 'bx')
        else
            plot(phi(i) - 2*pi, filtered_equilibria{i}(j, 1), 'rx')
        end 
        
    end

    figure(f4)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot(phi(i) - 2*pi, filtered_equilibria{i}(j, 2), 'bx')
        else
            plot(phi(i) - 2*pi, filtered_equilibria{i}(j, 2), 'rx')
        end 
        
    end
end


figure(f3)
grid on
xlabel("\theta_{r, d} [rad]")
ylabel("\theta_0 [rad]")
% legend("Unstable", "Stable")
xlim([-pi, pi])
title("Equilibria of Zero Dynamics: \theta_0")
hold on
xl1 = xline(-pi/2,'--',{'Blue-Sky', 'Catastrophe'}, 'Color', [0.4660 0.6740 0.1880]);
% xl1 = xline(-pi/2,'--',{'-\pi/2'}, 'Color', [0.4660 0.6740 0.1880]);
xl1.LabelVerticalAlignment = 'middle';
xl1.LabelHorizontalAlignment = 'center';

xl2 = xline(pi/2,'--',{'Blue-Sky', 'Catastrophe'}, 'Color', [0.4660 0.6740 0.1880]);
% xl2 = xline(pi/2,'--',{'\pi/2'}, 'Color', [0.4660 0.6740 0.1880]);
xl2.LabelVerticalAlignment = 'middle';
xl2.LabelHorizontalAlignment = 'center';
hold off

figure(f4)
grid on
xlabel("\theta_{r, d} [rad]")
ylabel("\theta_1 [rad]")
% legend("Unstable", "Stable")
xlim([-pi, pi])
title("Equilibria of Zero Dynamics: \theta_1")
hold on
xl1 = xline(-pi/2,'--',{'Blue-Sky', 'Catastrophe'}, 'Color', [0.4660 0.6740 0.1880]);
% xl1 = xline(-pi/2,'--',{'-\pi/2'}, 'Color', [0.4660 0.6740 0.1880]);
xl1.LabelVerticalAlignment = 'middle';
xl1.LabelHorizontalAlignment = 'center';

xl2 = xline(pi/2,'--',{'Blue-Sky', 'Catastrophe'}, 'Color', [0.4660 0.6740 0.1880]);
% xl2 = xline(pi/2,'--',{'\pi/2'}, 'Color', [0.4660 0.6740 0.1880]);
xl2.LabelVerticalAlignment = 'middle';
xl2.LabelHorizontalAlignment = 'center';
hold off

%% 3D Equilibria
f5 = figure;

for i = 1:length(phi)/2
    figure(f5)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot3(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), phi(i), 'bx')
        else
            plot3(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), phi(i), 'rx')
        end 
        
    end
end

for i = 32:length(phi)
    figure(f5)
    hold on
    for j = 1:size(filtered_equilibria{i}, 1)
        
        if(stability{i}(j))
            plot3(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), phi(i)-2*pi, 'bx')
        else
            plot3(filtered_equilibria{i}(j, 1), filtered_equilibria{i}(j, 2), phi(i)-2*pi, 'rx')
        end 
        
    end
end

% % % Add Plane
% [x y] = meshgrid(-15:0.1:15); % Generate x and y data
% z1 = pi/2 + zeros(size(x, 1)); % Generate z data
% s1 = surf(x, y, z1) % Plot the surface
% s1.EdgeColor = 'None';
% s1.FaceAlpha = 0.1;
% s1.FaceColor = [0.4660 0.6740 0.1880];
% 
% z2 = -pi/2 + zeros(size(x, 1)); % Generate z data
% s2 = surf(x, y, z2) % Plot the surface
% s2.EdgeColor = 'None';
% s2.FaceAlpha = 0.1;
% s2.FaceColor = [0.4660 0.6740 0.1880];
% % hold off

% load("result_step.mat");
load("result_PP.mat");
plot3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1))
quiver3(result.simout.data(:, 2), result.simout.data(:, 3), result.simout.data(:, 1), ...
        result.simout3.data(:, 2), result.simout3.data(:, 3), result.simout3.data(:, 1),'Color', [0.9290 0.6940 0.1250])
plot3(result.simout.data(1, 2), result.simout.data(1, 3), result.simout.data(1, 1), 'go')
hold off

figure(f5)
grid on
xlabel("\theta_0 [rad]")
ylabel("\theta_1 [rad]")
zlabel("\theta_{rd} [rad]")
legend("Unstable", "Stable")
% xlim([-pi, pi])
title("Phase Space with Equilibria: PP planner")

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