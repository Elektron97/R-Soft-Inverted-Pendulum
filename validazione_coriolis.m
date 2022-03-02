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
        
        testCr = coriolisMatrix(THETAR(i, j), theta0_sing, theta1_sing, 0, 0, 0, m, L, D);
        
        Cr11(i, j) = testCr(2, 2);
        Cr12(i, j) = testCr(2, 3);
        Cr21(i, j) = testCr(3, 2);
        Cr22(i, j) = testCr(3, 3);
        
    end
end

%% Plot

figure
subplot(2, 2, 1)
s1 = surf(THETA0, THETA1, real(Cr11));
s1.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{22}");
view(23, 21);

subplot(2, 2, 2)
s2 = surf(THETA0, THETA1, real(Cr12));
s2.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{23}");
view(23, 21);

subplot(2, 2, 3)
s3 = surf(THETA0, THETA1, real(Cr21));
s3.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{32}");
view(23, 21);

subplot(2, 2, 4)
s4 = surf(THETA0, THETA1, real(Cr22));
s4.EdgeColor = 'none';
xlabel("\theta_0");
ylabel("\theta_1");
zlabel("myBr_{33}");
view(23, 21);