function D = dampingMatrix(beta,beta_r)
%DAMPINGMATRIX
%    D = DAMPINGMATRIX(BETA,BETA_R)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    22-Feb-2022 10:17:39

t2 = beta./2.0;
D = reshape([beta_r,0.0,0.0,0.0,beta,t2,0.0,t2,beta./3.0],[3,3]);
