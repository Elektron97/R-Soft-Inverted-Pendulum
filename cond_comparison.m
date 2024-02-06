%%% Cond Number Comparison %%%
clear all
close all
clc

%% Load Datasets
load("condition_numberCP.mat")

cd ..
cd Master-Thesis/MATLAB/Renda_Model/RSoft_InvertedPendulum
load("cond_numberSP.mat")

cd ../../../../R-Soft-Inverted-Pendulum

%% Useful Vars
colormap_type = "turbo";

%% Flat Comparison
% % Condition Number of Inertia Matrix
% s5 = surfc(THETA0, THETA1, condBr_cp);
% s5(1).FaceColor = "#0072BD";
% % s5.LineStyle = ":";
% s5(1).EdgeColor = 'none';
% xlabel("\theta_0, q_0");
% ylabel("\theta_1, q_1");
% zlabel("log(\chi(M))");
% view(48, 13);
% % title("Condition Number of Inertia Matrix (CP)")
% set(gca, 'ZScale', 'log')
% % zlim([0, 1e+6]);
% % colorbar
% set(gca, 'ColorScale', 'log')
% hold on
% s6 = surfc(Q0, Q1, condBr);
% s6(1).EdgeColor = 'none';
% % s6.LineStyle = ":";
% s6(1).FaceColor = "#D95319";
% hold off
% % legend("Curv. Param.", "Strain Param.")

%% Curve Levels
% Condition Number of Inertia Matrix
figure
% Figure params
face_alpha = 0.8;
contour_width = 1.8;
contour_style = '-';
contours_cp = 1:15;
% contours_cp = 15;
contours_sp = [1 2 3 4 5 6];

% CP case
s5 = surf(THETA0, THETA1, log10(condBr_cp));
s5.EdgeColor = 'none';
s5.FaceColor = "#0072BD";
% s5.FaceColor = "#77AC30";
s5.FaceAlpha = face_alpha;
hold on
% SP case
s6 = surf(Q0, Q1, log10(condBr));
s6.EdgeColor = 'none';
s6.FaceColor = "#D95319";
% s6.FaceColor = "#D95319";
s6.FaceAlpha = face_alpha;

% Contours
contour3(THETA0, THETA1, log10(condBr_cp), contours_cp, 'LineWidth', contour_width, 'LineStyle', contour_style);
contour3(Q0, Q1, log10(condBr), contours_sp, 'LineWidth', contour_width, 'LineStyle', contour_style);
hold off
xlabel("\theta_0, q_0");
ylabel("\theta_1, q_1");
zlabel("log_{10}(\chi(M))");
colormap(colormap_type);
c = colorbar;
c.Label.String = 'Cond. Number \chi(M)';
c.Ticks = 2:3:14;
c.TickLabels = compose('10^{%d}',c.Ticks);
view(38, 5);
leg1 = legend("Curv. Param.", "Strain Param.");
set(leg1, 'Position',[0.131181547619048 0.717460317460319 0.223214285714286 0.0837301587301587]);

%% pcolor
% figure
% p1 = pcolor(THETA0, THETA1, log10(condBr_cp));
% p1.EdgeColor = 'none';
% shading interp;
% % set(gca,'ColorScale','log');
% colormap("turbo");
% colorbar
% 
% figure
% p2 = pcolor(Q0, Q1, log10(condBr));
% p2.EdgeColor = 'none';
% shading interp;
% colormap("turbo")
% colorbar

%% contourf
% figure
% pcolor(THETA0, THETA1, log10(condBr_cp));
% hold on
% shading interp;
% contour(THETA0, THETA1, log10(condBr_cp), [1 3 5 6 7 11 14], 'LineColor', 'k');
% colormap("turbo");
% colorbar

% figure
% contourf(THETA0, THETA1, log10(condBr_cp), [1 3 5 6 7 11 14]);
% colormap("turbo");
% colorbar
%% contour
% figure
% contour(THETA0, THETA1, log10(condBr_cp), 15);
% colormap("turbo");
% grid on
% colorbar

% figure
% pcolor(THETA0, THETA1, log10(condBr_cp));
% hold on
% shading interp;
% contour(THETA0, THETA1, log10(condBr_cp), [1 3 5 6 7 11 14], 'LineColor', 'k');
% colormap("turbo");
% colorbar