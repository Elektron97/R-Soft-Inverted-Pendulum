%%% Multiple Frame Plot %%%
clear all
close all
clc

%% Load Result
load("result_PP_with_obj.mat")

%% Useful variables
is_pp = true;
L = 1;
D = 0.1;

% Useful colors
gray_col = [0.7, 0.7, 0.7];
htmlGray = [128 128 128]/255;
purple_col = [0.4940 0.1840 0.5560];
orange_col = [0.8500 0.3250 0.0980];
green_col = [0.4660 0.6740 0.1880];
celeste_figo = [0 0.4470 0.7410];
cool_yellow = [0.9290 0.6940 0.1250];
cool_red = [0.6350 0.0780 0.1840];
black = [0 0 0];

%% Plot
% Insert skip of robot frame configuration
skip_frame_robot = 5;
% skip_frame_quiver = 7;

% Unique equilibrium
theta_r0 = 1.7;

theta_r1 = -(2/3)*pi - 0.1; % m_obj = 0.3
speed = -0.1;
t1 = 40;
t2 = (theta_r1 - theta_r0)/speed + t1;

idx1 = find(result.simout.time > t1);
idx1 = idx1(1) - 1;

idx2 = find(result.simout.time > t2);
idx2 = idx2(1) + 2;

if is_pp
    phase_time = [1, idx1, idx2, length(result.simout.time)];
    phase_color = {green_col, orange_col, cool_yellow, celeste_figo};

    %%%%% Cicle on Phases
    for phase=1:3

        figure
        hold on
        % Initial Position
        plot_Rsoft(result.simout.data(phase_time(phase), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase})
        % Final Position
        plot_Rsoft(result.simout.data(phase_time(phase+1), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase+1})
        for i = phase_time(phase):phase_time(phase+1)
            % skip_frame
            % rewrite (i % skip_frame_robot == 0) in matlab
            if(mod(i, skip_frame_robot) == 0)
                plot_Rsoft(result.simout.data(i, :), L, D, plot_frame=false, plot_thick=false, color=gray_col, backbone_thick=1)
            end
            % Trajectory always continous
            tip_traj(:, i) = fwdKinRSIP(result.simout.data(i, :), 1, 0, L, D);
            % tip_vel(:, i) = jacobianSoft(result.simout.data(i, 1), ...
            %                                 result.simout.data(i, 2), ...
            %                                 result.simout.data(i, 3), ...
            %                                 1, 0, L, D)*result.simout3.data(i, :)';
        end
        % Initial Position
        plot_Rsoft(result.simout.data(phase_time(phase), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase})
        % Final Position
        plot_Rsoft(result.simout.data(phase_time(phase+1), :), L, D, plot_frame=false, plot_thick=false, color=phase_color{phase+1})
        hold on
        % Tip Trajectory
        plot3(tip_traj(1, phase_time(phase):phase_time(phase+1)), ...
            tip_traj(2, phase_time(phase):phase_time(phase+1)), ...
            zeros(1, length(phase_time(phase):phase_time(phase+1))), ...
            'LineWidth', 2, 'Color', htmlGray)
        %Tip Velocity: all points
        % quiver3(tip_traj(1, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
        %         tip_traj(2, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
        %         zeros(1, length(phase_time(phase):skip_frame_quiver:phase_time(phase+1))), ...
        %         tip_vel(1, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
        %         tip_vel(2, phase_time(phase):skip_frame_quiver:phase_time(phase+1)), ...
        %         zeros(1, length(phase_time(phase):skip_frame_quiver:phase_time(phase+1))), ...
        %         'Color', '#25283D', 'LineWidth', 1.0)

        % Add Object
        % TO DO
        
        hold off
        xlabel("x [m]")
        ylabel("y [m]")
        % title("Pick and Place: Phase " + num2str(phase))
        legend("Initial Config. of Phase " + num2str(phase), "Final Config. of Phase " + num2str(phase))

        clear tip_traj
        clear tip_vel
    end
else
    figure
    hold on
    for i =1:length(result.simout.time)
        plot_Rsoft(result.simout.data(i, :), L, D, plot_frame=false, plot_thick=false, color=gray_col, backbone_thick=1)
        tip_traj(:, i) = fwdKinRSIP(result.simout.data(i, :), 1, 0, L, D);
    end
    % Initial Position
    plot_Rsoft(result.simout.data(1, :), L, D, plot_frame=false, plot_thick=false, color=green_col)
    % Final Position
    plot_Rsoft(result.simout.data(end, :), L, D, plot_frame=false, plot_thick=false, color=orange_col)
    hold on
    % Tip Trajectory
    plot3(tip_traj(1, :), tip_traj(2, :), zeros(1, length(result.simout.time)), ...
        'LineWidth', 2, 'Color', htmlGray)
    hold off
    xlabel("x [m]")
    ylabel("y [m]")
end