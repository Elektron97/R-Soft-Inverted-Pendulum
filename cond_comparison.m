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
% % Condition Number of Inertia Matrix
% figure
% % Figure params
% face_alpha = 0.7;
% contour_width = 1.5;
% contour_style = '-';
% 
% % CP case
% s5 = surf(THETA0, THETA1, log10(condBr_cp));
% s5.EdgeColor = 'none';
% s5.FaceColor = "#0072BD";
% s5.FaceAlpha = face_alpha;
% hold on
% % SP case
% s6 = surf(Q0, Q1, log10(condBr));
% s6.EdgeColor = 'none';
% s6.FaceColor = "#D95319";
% s6.FaceAlpha = face_alpha;
% 
% % Contours
% contour3(THETA0, THETA1, log10(condBr_cp), 15, 'LineWidth', contour_width, 'LineStyle', contour_style);
% contour3(Q0, Q1, log10(condBr), 5, 'LineWidth', contour_width, 'LineStyle', contour_style);
% hold off
% xlabel("\theta_0, q_0");
% ylabel("\theta_1, q_1");
% zlabel("log_{10}(\chi(M))");
% colormap("turbo");
% c = colorbar;
% c.Label.String = 'Cond. Number \chi(M)';
% c.Ticks = 2:3:14;
% c.TickLabels = compose('10^{%d}',c.Ticks);
% view(48, 13);
% leg1 = legend("Curv. Param.", "Strain Param.");
% set(leg1, 'Position',[0.131181547619048 0.717460317460319 0.223214285714286 0.0837301587301587]);

%% pcolor
% figure
% p1 = pcolor(THETA0, THETA1, log10(condBr_cp));
% p1.EdgeColor = 'none';
% % set(gca,'ColorScale','log');
% colorbar
% 
% figure
% p2 = pcolor(Q0, Q1, log10(condBr));
% p2.EdgeColor = 'none';
% colorbar

%% contourf
% figure
% contourf(THETA0, THETA1, log10(condBr_cp), [1, 2, 5, 6, 9, 14]);
% colorbar