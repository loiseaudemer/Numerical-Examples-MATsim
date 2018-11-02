clear all
clc


A = [0.8 0;0 -1.2];
B = [0;1];
Q = [2 0;0 1];
R = 1;

nu = 1;
nx = 2;

% Laplacian matrix
% LD = diag([1,2,2,2,2,2,2,1]);
% LA = [0 1 0 0 0 0 0 0;
%       1 0 1 0 0 0 0 0;
%       0 1 0 1 0 0 0 0;
%       0 0 1 0 1 0 0 0;
%       0 0 0 1 0 1 0 0;
%       0 0 0 0 1 0 1 0;
%       0 0 0 0 0 1 0 1;
%       0 0 0 0 0 0 1 0];
%   
% L = LD - LA;

LD4 = diag([1,2,2,1]);
LA4 = [0 1 0 0;
      1 0 1 0;
      0 1 0 1;
      0 0 1 0];

  
  

L = LD4 - LA4;


% eigenvalues of L
lambdaM = max(eig(L));
lambda2 = sort(eig(L));
lambda2 = lambda2(2);

c = 2/(lambda2 + lambdaM);
cl = -(c^2*lambdaM^2 - 2*c*lambdaM);
Rcl = sqrt(cl*R^(-1));
Bcare = B*Rcl;
% Find K and Q
epsilon = 1e-3;
Qcare = lambdaM*Q + epsilon*eye(2);

[P,PL,PG] = dare(A,B,Qcare,Rcl);
%K = -c*inv(R+B'*P*B)*B'*P*A;
K = -c*inv(R+B'*P*B)*B'*P*A;

Phi = A + B*K;


%% Initial conditions
agent = 4;


x10 = [-0.08;0.11];
x20 = [0.12;-0.08];
x30 = [-0.09;-0.14];
x40 = [-0.12;0.04];
% x50 = [0.07;-0.16];
% x60 = [-0.21;0.12];
% x70 = [0.15;-0.22];
% x80 = [-0.17;-0.14];

x0 = [x10;x20;x30;x40];



%% Distribute Model Predictive Control Method
mpciter = 2;
N = 10; % Prediction horizon



Xtraj = x0;




for iter = 1:mpciter
    
    % state estimation for initial step
    if iter == 1    % first step estimation using x+=Phi*x
        Xestimate(:,1) = x0;
        for  agentCount = 1:agent
            for pred_iter = 1:N %(0:N-1)
                Xestimate([(2*agentCount-1):(2*agentCount)],pred_iter+1) = ...
                    Phi*Xestimate([(2*agentCount-1):(2*agentCount)],pred_iter);
                    %Xestimate([(2*agentCount-1):(2*agentCount)],pred_iter+1) = zeros(2,1);
            end
        end
    end
    
    
    
    for agentOpt = 1:agent
        
        ui = sdpvar(repmat(nu,1,N),repmat(1,1,N));
        xi = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));
        ci = sdpvar(repmat(nu,1,N),repmat(1,1,N));
        
        
        costi = 0;
        constraints = [xi{1} == Xtraj([(2*agentOpt-1):(2*agentOpt)],iter)];

        
        % sum cost(0:N-1) and constraints(0:N-1)
        for n = 1:N

            % find \sum(xi - xj)
            if agentOpt == 1 % x1 - x2
                Xdisa = xi{n} - Xestimate([3:4],n);
            elseif agentOpt == 2 % x2 -x1 + x2 - x3
                Xdisa = 2*xi{n} - Xestimate([1:2],n) - Xestimate([5:6],n);
            elseif agentOpt == 3 % 2*x3 - x2 - x4
                Xdisa = 2*xi{n} - Xestimate([3:4],n) - Xestimate([7:8],n);
            elseif agentOpt == 4
                Xdisa = xi{n} - Xestimate([5:6],n); 
            end
            
            
            
            ui{n} = ci{n} + K*Xdisa;
            costi = costi +  ci{n}'*2*ci{n}; %+ Xdisa'*eye(2)*Xdisa; %ui{n}'*2*ui{n}; %Xdisa'*eye(2)*Xdisa;
            xi{n+1} = A*xi{n} + B*ui{n};
        end
        
        constraints = [constraints, ci{N} == 0];
        
        ops = sdpsettings('solver','mosek');
        optimize(constraints,costi,ops);
        
        
        % implement the opt input to agent i
        ui_star = value(ui{1});
        
        % record real systems trajs
        Utraj(agentOpt,iter) = ui_star;
        UKXtraj(agentOpt,iter) = K*Xtraj([2*agentOpt,2*agentOpt-1],iter);
        Xtraj([(2*agentOpt-1):(2*agentOpt)],iter+1) = ...
            A*Xtraj([(2*agentOpt-1):(2*agentOpt)],iter) + B*ui_star;
        
        %record the optimal state sequence and input sequence, Estimation
        %for next iter
        for i = 1:N
            XESTIMATE([(2*agentOpt-1):(2*agentOpt)],i) = value(xi{i});
        end
        
    end
        XESTIMATE(:,N+1) = kron(eye(4),Phi)*XESTIMATE(:,N);
    % update Xestimate
    XESTIMATE = XESTIMATE(:,[2:end]);
    Xestimate = XESTIMATE;
        
end







% Test with systems

% Test for the Integrated System
% Xtraj = x0;
% 
% Xestimate = x0;
% 
% 
% 
% for k = 1:100
%     
%     for agentComp = 1:agent
%         
%         if (agentComp == 1)
%             Xdisa = Xtraj([2*agentComp-1,2*agentComp],k) - Xestimate([3:4],1);
%         elseif agentComp == 2
%             Xdisa = 2*Xtraj([2*agentComp-1,2*agentComp],k) - Xestimate([1:2],1) - Xestimate([5:6],1);
%         elseif agentComp == 3
%             Xdisa = 2*Xtraj([2*agentComp-1,2*agentComp],k) - Xestimate([3:4],1) - Xestimate([7:8],1);
%         elseif agentComp == 4
%             Xdisa = Xtraj([2*agentComp-1,2*agentComp],k) - Xestimate([5:6],1); 
%         end
%         
%         
%         Utraj(agentComp,k) = K*Xdisa;
%         Xtraj([(2*agentComp-1):2*agentComp],k+1) = A*Xtraj([(2*agentComp-1):2*agentComp],k) +...
%             B*Utraj(agentComp,k);
%     end
%         Xestimate = Xtraj(:,k);
%         
%  
% end


% Test for single system


% % 
% Xtraj = x0;
% for k = 1:mpciter
%     
%     Utraj(:,k) = kron(L,K)*Xtraj(:,k);
%     Xtraj(:,k+1) = kron(eye(agent),A)*Xtraj(:,k) + kron(eye(agent),B)*Utraj(:,k);
%     
% end

%% Plot trajectories



steps = 1:mpciter+1;


figure(1);
for i = 1:2:(2*agent)
    scatter(steps, Xtraj(i,:));
    hold on 
    plot(steps, Xtraj(i,:),'-.');
    hold on
end

ylabel('x_1')
xlabel('steps')
grid on 




figure(2);
for i = 2:2:(2*agent)
    scatter(steps, Xtraj(i,:));
    hold on
    plot(Xtraj(i,:),'-.');
    hold on
end
ylabel('x_2')
xlabel('steps')
grid on


figure(3);
for i = 1:2:2*agent
    plot(Xtraj(i,:),Xtraj(i+1,:),'-.');
    hold on
end
xlabel('x_1');
ylabel('x_2');
grid on


figure(4);
for i = 1:agent
    stairs(Utraj(i,:));
    hold on
end
xlabel('steps')
ylabel('u')
grid on


% figure(4);
% for i = 1:agent
%     stairs(UKXtraj(i,:));
%     hold on
% end
% xlabel('steps')
% ylabel('u')
% grid on
