%% Clear workspace
yalmip('clear')
clear all
clc

%% Dynamics

A = [1.1 1;0 1.3];
B = [1;1];
D = eye(2);     % disturbance mat
nu = 1;
nx = 2;
nr = 1;
%% linear LQ feedback 

Q = eye(2,2);
R = 0.01;

[F,S,e] = dlqr(A,B,Q,R,0);
F = -F;

%% closed-loop system parameters, x+ = Phi*x + D*w

Phi = A + B*F;

%% find \Psi by solving discrete-time riccati equation

[P,L,G] = dare(A,B,Q,R);
Psi = R + B'*P*B;


%% Prediction

% mpc, nominal
N = 10;
x0 = [-6.9;2.38];
%x0 = [1;-2];
mpciter = 20;
ops = sdpsettings('solver','mosek');




% one step

for iter = 1:mpciter
    
    x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
    c = sdpvar(repmat(nr,1,N),repmat(1,1,N));
    
    constraints = [];
    obj = 0;
    
    constraints = [x{1} == x0];
    for k = 1:N
        obj = obj + c{k}'*Psi*c{k};
        
        constraints = [constraints, -1 <= F*x{k} + c{k} <= 1];
        constraints = [constraints, -10 <= x{k} <= 10];
        x{k+1} = Phi*x{k} + B*c{k};
        
    end
    constraints = [constraints, -10 <= x{N+1} <= 10];
    
    optimize(constraints,obj,ops);
    
    
    controller(:,iter) = value(F*x{1} + c{1});
    
    dist = (rand(2,1) - 0.5*ones(2,1)) * (0.12/0.5);
    %dist = [-0.12;-0.12];
    x0 = A*x0 + B*controller(:,iter) + D*dist;
    Xtraj(:,iter) = x0;
    
end

Xtraj = [[-6.9;2.38] Xtraj];


figure(1);
scatter(Xtraj(1,:),Xtraj(2,:));
hold on 
plot(Xtraj(1,:),Xtraj(2,:));
hold on
plot(plotRSet(Phi,D,0.12,20),[],'red');
%axis([-8 8 -4 4]);
hold off
grid on

figure(2);
stairs(controller);


