%% This code is to find control singnals when l2norm of the terminal state are 0.0870 in conventional method and proposed method respectively
clear

%% System model(rocket)
% Plant matrixes
A = [0 1; 0 0];
b = [0; 1];
% System size (Number of state)
d = length(b);
% Initial state
x0 = [10; -3];
% Horizon lengths(terminal time)
T = 5.0;
%% Time discretization
% Discretization size
n = 100; 
% Discretization interval
h = T/n; 
% System discretization
[Ad,bd] = c2d(A,b,h);
% Matrix Phi
Phi = zeros(d,n);
v = bd;
Phi(:,end) = v;
for j = 1:n-1
    v = Ad*v;
    Phi(:,end - j) = v;
end
% Vector zeta
zeta = -Ad^n*x0;
% Phi*u-zeta is the vector discretize terminal state

%% Simulation preparation 
% Number of simulation
MAX_ITER = 5000;
% l2 norm of terminal states are 0.0870
% Hyperparameters othar than tilde lambda 
gamma = 0.2;
alpha = 1;
beta = 1; 
% Lambda (regularization parameter) which is used in conventional method
tilde_lambda_con = 0.0087;
% Lambda which is used in proposed method
tilde_lambda_pro = 0.0218;
% Matrix Psi
Psi = [Phi; eye(n); eye(n)];
% Matrix M (matrix which is used to find u)
M = (Psi'*Psi) \ Psi';
% Dimension number of Psi
m = d + 2*n;
% Saturation function
sat = @(x) sign(x).*min(abs(x),1);

% Initial vectors in conventional method
% u1 is control signal which is generated in conventional method
u1 = zeros(n, 1); z1 = zeros(m, 1); v1 = zeros(m, 1); 
% Initial vectors in proposed method
% u2 is control signal which is generated in proposed method
u2 = zeros(n, 1); z2 = zeros(m, 1); v2 = zeros(m, 1); 

% Value which is used to find z
gamma_term = 1/(2*gamma + 1);
%% Simulation
tic
% Repeat MAX_ITER times
for k = 1:MAX_ITER
    % Update u,z,v in conventional method
    u1 = M*(z1 - v1);       
    z1(1:d) = gamma_term*(2*gamma* zeta + (Phi*u1 + v1(1:d)));
    z1(d+1:n+d) = soft_thresholding(u1 + v1(d+1:n+d), gamma*tilde_lambda_con);
    z1(n+d+1:m) = sat(u1 + v1(n+d+1:m));
    v1 = v1 + Psi*u1 - z1;
end
% Computation time in conventional method
com_con = toc;

tic
% Repeat MAX_ITER times
for k = 1:MAX_ITER
    % Update u,z,v in proposed method
    u2 = M*(z2 - v2);
    z2(1:d) = gamma_term*(2*gamma*zeta + (Phi*u2 + v2(1:d)));
    z2(d+1:n+d) = firm_thresholding(u2 + v2(d+1:n+d), alpha*gamma*tilde_lambda_pro, alpha*beta);
    z2(n+d+1:m) = sat(u2 + v2(n+d+1:m));
    v2 = v2 + Psi*u2 - z2;
end
% Computation time in proposed method
com_pro = toc;

% Terminal state in conventional method
terminal_state_con = Phi*u1 - repmat(zeta, 1, 1);
% Terminal state in proposed method
terminal_state_pro = Phi*u2 - repmat(zeta, 1, 1);

% l2 norm of terminal state in conventional method
l2norm_con = sqrt(sum(terminal_state_con.^2, 1))';
% l2 norm of terminal state in proposed method
l2norm_pro = sqrt(sum(terminal_state_pro.^2, 1))';

% Regarded as non-zero when the control signal is 0.00001 or greater
% l0 norm of u (control signal) in conventional method 
l0norm_con = sum(abs(u1) >= 0.00001, 1)';
% l0 norm of u (control signal) in proposed method 
l0norm_pro = sum(abs(u2) >= 0.00001, 1)';

% Sparsity (Number of discrete time intervals where the control signal is 0)
% Sparsity in conventional method 
sparsity_con = n - l0norm_con;
% Sparstiy in proposed method
sparsity_pro = n - l0norm_pro;

% Display each computaion time
com_con_time = sprintf('Computation time in conventional method: %.4f (s)', com_con);
com_pro_time = sprintf('Computation time in proposed method: %.4f (s)', com_pro);
disp(com_con_time);
disp(com_pro_time);

% Display each l2 norm 
l2norm_con_val = sprintf('l2 norm of terminal state in conventional method: %.4f', l2norm_con);
l2norm_pro_val = sprintf('l2 norm of terminal state in proposed method: %.4f', l2norm_pro);
disp(l2norm_con_val);
disp(l2norm_pro_val);

% Display each sparsity
sparsity_con_val = sprintf('Sparsity in conventional method: %d', sparsity_con);
sparsity_pro_val = sprintf('Sparsity in proposed method: %d', sparsity_pro);
disp(sparsity_con_val);
disp(sparsity_pro_val);

%% plot control signals
figure;
% plot control signal in conventional method
plot(0:T/n:T-T/n,u1,'--','LineWidth',3);
hold on 
% plot control signal in proposed method
plot(0:T/n:T-T/n,u2,'LineWidth',3);
hold off
grid on
grid minor
legend('l1 norm','mcpf','NumColumns',2);
xlabel('Times (s)');
ylabel('Control signals');
title('The control signals when l2norm of terminal state is 0.087');
