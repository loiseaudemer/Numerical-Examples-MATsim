%% Numerical example simulation
%   A suboptimal approach to distributed linear quadratic

clear all
clc

% dynamics of the agents
agent = 8;
A = [0 1;-1 0];
B = [0;1];

% number of states and inputs
nx = 2;
nu = 1;

% Laplacian matrix
LD = diag([1,2,2,2,2,2,2,1]);
LA = [0 1 0 0 0 0 0 0;
      1 0 1 0 0 0 0 0;
      0 1 0 1 0 0 0 0;
      0 0 1 0 1 0 0 0;
      0 0 0 1 0 1 0 0;
      0 0 0 0 1 0 1 0;
      0 0 0 0 0 1 0 1;
      0 0 0 0 0 0 1 0];
  
L = LD - LA;

% eigenvalues of L
lambda8 = max(eig(L));
lambda2 = sort(eig(L));
lambda2 = lambda2(2);

% cost function
Q = [2 0;0 1];
R = 1;


% suboptimal bound \gamma
gamma = 3;

% find P by solving the Riccati equation
epsilon = 1e-3;
c = 2/(lambda2 + lambda8);
cl = -(c^2*lambda8^2 - 2*c*lambda8);
Rcl = sqrt(cl*R^(-1));
Bcare = B*Rcl;
Qcare = lambda8*Q + epsilon*eye(2);

[P,PL,PG] = care(A,Bcare,Qcare);

% control gain K
K = -c*R^(-1)*B'*P;


%% Initial conditions
x10 = [-0.08;0.11];
x20 = [0.12;-0.08];
x30 = [-0.09;-0.14];
x40 = [-0.12;0.04];
x50 = [0.07;-0.16];
x60 = [-0.21;0.12];
x70 = [0.15;-0.22];
x80 = [-0.17;-0.14];

x0 = [x10;x20;x30;x40;x50;x60;x70;x80];


%% Iters with the decoupled control law
AA = diag(sort(svd(L)));
[U,eigval] = eig(L);

% integrated system dynamics
Abar = (kron(eye(agent),A) + kron(AA,B*K));

x0b = kron(U',eye(nx))*x0;


T = 0.01;
xb = x0b;
for k = 1:1500
    xb(:,k+1) = (eye(16) + T*Abar)*xb(:,k);
end

x = inv(kron(U',eye(nx)))*xb;




% plot x_1
figure(1);
for traj = 1:2:2*agent
    plot(x(traj,:));
    hold on 
end
hold off
grid on
xlabel('steps')
ylabel('x_1')
axis([0 1700 -0.25 0.17]);





% plot x_2
figure(2);
for traj = 2:2:2*agent
    plot(x(traj,:));
    hold on 
end
hold off
grid on
axis([0 1700 -0.25 0.2]);
xlabel('steps')
ylabel('x_2')

% % plot trajs
% figure(3);
% for traj = 1:agent
%     plot(x(traj,:),x(traj+1,:));
%     hold on
%     grid on
% end
% xlabel('x_1');
% ylabel('x_2');


