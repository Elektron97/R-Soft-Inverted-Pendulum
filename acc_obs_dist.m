%%%%%%%%%%%% Accessibility and Observation Distribution %%%%%%%%%%%%
%% Compute Model
Rsoft_model

%% Accessibility Distribution
% LfG = lieBracket(F, G, x);
% disp("First Order Computed");
% L2fG = lieBracket(F, LfG, x);
% LgLfG = lieBracket(G, LfG, x);
% disp("Second Order Computed");
% L3fG = lieBracket(F, L2fG, x);
% LfLgLfG = lieBracket(F, LgLfG, x);
% % Skip other lieBracket
% disp("Third Order Computed");

%% Observability Codistribution
% Possible outputs
h1 = theta_r;
dh = jacobian(h1, x);

disp("Observability codistribution with outputs: ")
h1

% First Order
Lfdh = rowLieBracket(dh, F, x);
% Lgdh = rowLieBracket(dh, G, x); % Zeros
disp("First order computed");

% Second Order
Lf2dh = rowLieBracket(Lfdh, F, x);
% Lg2dh = rowLieBracket(Lgdh, G, x);
LgLfdh = rowLieBracket(Lfdh, G, x);
% LfLgdh = rowLieBracket(Lgdh, F, x);
disp("Second order computed");

% Third Order
Lf3dh = rowLieBracket(Lf2dh, F, x);
Lg2Lfdh = rowLieBracket(LgLfdh, G, x);
disp("Third Order computed");

obs_codistr = [dh; Lfdh; Lf2dh; LgLfdh; Lf3dh; Lg2Lfdh];
disp("Save codistribution");
save('obs_codistr.mat', 'obs_codistr');

% % Possible outputs
% h1 = theta_r;
% dh = jacobian(h1, x);
% 
% disp("Observability codistribution with outputs: ")
% h1
% 
% % First Order
% Lfdh = rowLieBracket(dh, F, x);
% Lgdh = rowLieBracket(dh, G, x);
% disp("First order computed");
% 
% % Second Order
% Lf2dh = rowLieBracket(Lfdh, F, x);
% Lg2dh = rowLieBracket(Lgdh, G, x);
% LgLfdh = rowLieBracket(Lfdh, G, x);
% LfLgdh = rowLieBracket(Lgdh, F, x);
% disp("Second order computed");