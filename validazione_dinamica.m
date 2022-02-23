%%%%%%%%% Validazione Matrici Dinamiche %%%%%%%%%%%%%%
clear all
close all
clc

%% Add Functions
addpath("my_functions");
addpath(genpath("Della Santina"));
addpath("origin_soft_pendulum");
%% Parameters
L = 1;
D = 0.1;

m = 1;
g = 9.81;
beta = 0.1;
k = 1;

threshold = 1e-3;

%% Validation Gravity Vector
% [THETA0, THETA1] = meshgrid(-10:0.1:10, -10:0.1:10);
% THETAR = 0*THETA0;
% 
% for i = 1:size(THETA0, 1)
%     for j = 1:size(THETA0, 2)
%         
%         % Not Singular Configurations
%         if(abs(THETA1(i, j)) <= threshold)
%            theta1_sing = threshold;
%            
%            if(abs(THETA0(i, j)) <= threshold)
%                 theta0_sing = threshold;
%            else
%                 theta0_sing = THETA0(i, j);  
%            end   
%            
%         else
%             theta0_sing = THETA0(i, j);
%             theta1_sing = THETA1(i, j);
%         end
%         
%         testG = originGravityVector(theta0_sing, theta1_sing, m, g, L, D);
%         santinaG = G_fcn([theta0_sing; theta1_sing]);
%         testGr = gravityVector(THETAR(i, j), theta0_sing, theta1_sing, m, g, L, D);
%         
% %         testG = originGravityVector(THETA0(i, j), THETA1(i, j), m, g, L, D);
% %         santinaG = G_fcn([THETA0(i, j); THETA1(i, j)]);
% %         testGr = gravityVector(THETAR(i, j), THETA0(i, j), THETA1(i, j), m, g, L, D);
%         
%         G1(i, j) = testG(1);
%         G2(i, j) = testG(2);
%         
%         santina1(i, j) = santinaG(1);
%         santina2(i, j) = santinaG(2);
%         
%         Gr2(i, j) = testGr(2);
%         Gr3(i, j) = testGr(3);
%     end
% end
% 
% %G3
% figure
% s1 = surf(THETA0, THETA1, real(G1));
% s1.EdgeColor = 'none';
% xlabel("\theta_0");
% ylabel("\theta_1");
% zlabel("myG_1");
% 
% figure
% s2 = surf(THETA0, THETA1, real(santina1));
% s2.EdgeColor = 'none';
% xlabel("\theta_0");
% ylabel("\theta_1");
% zlabel("Santina G1");
% 
% figure
% s3 = surf(THETA0, THETA1, real(Gr2));
% s3.EdgeColor = 'none';
% xlabel("\theta_0");
% ylabel("\theta_1");
% zlabel("myGr_1");
% 
% %G2
% 
% figure
% s1 = surf(THETA0, THETA1, real(G2));
% s1.EdgeColor = 'none';
% xlabel("\theta_0");
% ylabel("\theta_1");
% zlabel("myG_2");
% 
% figure
% s2 = surf(THETA0, THETA1, real(santina2));
% s2.EdgeColor = 'none';
% xlabel("\theta_0");
% ylabel("\theta_1");
% zlabel("Santina G2");
% 
% figure
% s3 = surf(THETA0, THETA1, real(Gr3));
% s3.EdgeColor = 'none';
% xlabel("\theta_0");
% ylabel("\theta_1");
% zlabel("myGr_2");

%% Validation Inertia Matrix
[THETA0, THETA1] = meshgrid(-10:0.1:10, -10:0.1:10);
THETAR = 0*THETA0;


for i = 1:size(THETA0, 1)
    for j = 1:size(THETA0, 2)     
        
        % Not Singular Configurations
        if(abs(THETA1(i, j)) <= threshold)
           theta1_sing = threshold;
           
           if(abs(THETA0(i, j)) <= threshold)
                theta0_sing = threshold;
           else
                theta0_sing = THETA0(i, j);  
           end   
           
        else
            theta0_sing = THETA0(i, j);
            theta1_sing = THETA1(i, j);
        end
        
        
        testB = originInertiaMatrix(theta0_sing, theta1_sing, m, L, D);
        santinaB = B_fcn([theta0_sing; theta1_sing]);
        testBr = inertiaMatrix(THETAR(i, j), theta0_sing, theta1_sing, m, L, D);
        
%         testB = originInertiaMatrix(THETA0(i, j), THETA1(i, j), m, L, D);
%         santinaB = B_fcn([THETA0(i, j); THETA1(i, j)]);
%         testBr = inertiaMatrix(THETAR(i, j), THETA0(i, j), THETA1(i, j), m, L, D);
        
        B11(i, j) = testB(1, 1);
        B12(i, j) = testB(1, 2);
        B21(i, j) = testB(2, 1);
        B22(i, j) = testB(2, 2);
        
        santina11(i, j) = santinaB(1, 1);
        santina12(i, j) = santinaB(1, 2);
        santina21(i, j) = santinaB(2, 1);
        santina22(i, j) = santinaB(2, 2);
        
        Br11(i, j) = testBr(2, 2);
        Br12(i, j) = testBr(2, 3);
        Br21(i, j) = testBr(3, 2);
        Br22(i, j) = testBr(3, 3);
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% my B
figure
subplot(2, 2, 1)
s1 = surf(THETA0, THETA1, real(B11));
s1.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myB_{11}");
view(23, 21);

subplot(2, 2, 2)
s2 = surf(THETA0, THETA1, real(B12));
s2.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myB_{12}");
view(23, 21);

subplot(2, 2, 3)
s3 = surf(THETA0, THETA1, real(B21));
s3.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myB_{21}");
view(23, 21);

subplot(2, 2, 4)
s4 = surf(THETA0, THETA1, real(B22));
s4.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myB_{22}");
view(23, 21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Santina B
figure
subplot(2, 2, 1)
s1 = surf(THETA0, THETA1, real(santina11));
s1.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("Santina B_{11}");
view(23, 21);

subplot(2, 2, 2)
s2 = surf(THETA0, THETA1, real(santina12));
s2.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("Santina B_{12}");
view(23, 21);

subplot(2, 2, 3)
s3 = surf(THETA0, THETA1, real(santina21));
s3.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("Santina B_{21}");
view(23, 21);

subplot(2, 2, 4)
s4 = surf(THETA0, THETA1, real(santina22));
s4.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("Santina B_{22}");
view(23, 21);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% my Br
figure
subplot(2, 2, 1)
s1 = surf(THETA0, THETA1, real(Br11));
s1.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{22}");
view(23, 21);

subplot(2, 2, 2)
s2 = surf(THETA0, THETA1, real(Br12));
s2.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{23}");
view(23, 21);

subplot(2, 2, 3)
s3 = surf(THETA0, THETA1, real(Br21));
s3.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{32}");
view(23, 21);

subplot(2, 2, 4)
s4 = surf(THETA0, THETA1, real(Br22));
s4.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{33}");
view(23, 21);

