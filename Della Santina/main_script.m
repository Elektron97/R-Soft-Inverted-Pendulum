init_script

global p_vect

H = [1, 1/2; 1/2, 1/3];
iH = inv(H);


% Put initial condition here
% x_0 = [pi/4; -pi/4];


% evaluate control gains
p_vect = lqr(A_lin_fcn(m,g,L,k,beta),[0;0;1;1],diag([1;1;1;1]),1)

% call simulink
out = sim('dynamics_controlled');


%
% from here you have just plots
%

figure
plot(out.x_ev.Time,out.x_ev.Data, 'linewidth', 2)
h = legend({'$\theta_0$','$\theta_1$'});
set(h,'Interpreter','latex')
set(gca,'FontSize',18)
xlabel('Time $t$','interpreter','latex','FontSize',22)
ylabel('$\theta$ (rad)','interpreter','latex','FontSize',22)
set(gca,'TickLabelInterpreter','latex')
grid on, box on

figure
plot(out.x_ev.Time,out.tau_ev.Data, 'linewidth', 2)
set(gca,'FontSize',18)
xlabel('Time $t$','interpreter','latex','FontSize',22)
ylabel('$\tau$','interpreter','latex','FontSize',22)
set(gca,'TickLabelInterpreter','latex')
grid on, box on

figure
plot(out.x_ev.Time,(H*(out.x_ev.Data'))', 'linewidth', 2)
h = legend({'$\alpha(1)$','$\int_0^1 \alpha(s) \mathrm{d}s$'});
set(h,'Interpreter','latex')
set(gca,'FontSize',18)
xlabel('Time $t$','interpreter','latex','FontSize',22)
ylabel('$q$ (rad)','interpreter','latex','FontSize',22)
set(gca,'TickLabelInterpreter','latex')
grid on, box on

figure, plot_robot(out.x_ev)