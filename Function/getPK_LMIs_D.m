% Discrete-time LPV system
function [P,K] = getPK_LMIs_D(Q,R,constants)
% System dynamics
options = optimset('Display','on');
con =[];

% System dimensions
n = length(Q(1,:)); % state dimension
m = length(R(1,:)); % input dimension

% Weighting matrices for cost function
Q_eps = sqrtm(Q);%could include +epsilon here for nonlinear systems
R_eps = sqrt(R);
 
% X Y and Lambda matrices generation
Y_0_LMIs = sdpvar(m,n); 
x_i = sdpvar(n,1);
%p_i = sdpvar(n,1);
% Build Pbarxj
P_xj={};
for i = 1:constants.nx
    P_xj{1,i}=sdpvar(n);%sdpvar(n,n,'diagonal');
end
% Build Qbarxij
Q_xij={};
for p = 1:constants.pbar-1
    for i = 1:constants.nx
        Q_xij{p,i}=sdpvar(n);%,n,'diagonal');
    end
end
%% Multi Convexity constraints
% Starting to build the contraints
i = 0;
% Let's start by implementing (29)
% We check conditions only for the vertices of Theta at pbar of different
% states
for num_state = 1:constants.nx
    Theta_V = constants.FPS{num_state,constants.pbar}.V';
    for j=1:size(Theta_V,2)
        % Parameter vector p
          p = Theta_V(:,j);
          [A,B,C,D ] = getAiBiCiDi(p,constants,constants.pbar);
          Y = Y_0_LMIs; 
          X = diag(x_i);
          ineq = [x_i(num_state), (A*X+B*Y);
               (A*X+B*Y)', P_xj{1,num_state}];
          con = [con;ineq>=0];
          % Now 'con' is complete and we can start the optimization problem
           i=i+1;
    end
end
% Implementation of (30)
% We check conditions only for the vertices of Theta at different p and different
% states

for pred = 1:constants.pbar-1
    for num_state = 1:constants.nx
        Theta_V = constants.FPS{num_state,pred}.V';
        for j=1:size(Theta_V,2)
            % Parameter vector p
              p = Theta_V(:,j);
              [A,B,C,D] = getAiBiCiDi(p,constants,pred);
              Y = Y_0_LMIs; 
              X = diag(x_i);
              ineq = [1/Q(num_state,num_state), (C*X+D*Y);
                   (C*X+D*Y)', Q_xij{pred,num_state}];
              con = [con;ineq>=0];
              % Now 'con' is complete and we can start the optimization problem
               i=i+1;
        end
    end
end

% implementing (33)
PbarQ = zeros(n,n);
for state = 1:constants.nx
    for p = 1:constants.pbar-1
        PbarQ = PbarQ+Q_xij{p,state};
    end
    PbarQ = PbarQ+P_xj{1,state};
end

%implementing (32)
ineq = [X, PbarQ, (R_eps*Y)', (Q_eps*X)';...
       PbarQ, PbarQ, zeros(n,n+m);... 
       R_eps*Y, zeros(m,n), eye(m), zeros(m,n);...
      Q_eps*X, zeros(n,n+m), eye(n)];

con = [con;ineq>=0];


% Display the number of iterations
disp('Number of iterations to get constraints for optimization to get P and K');
disp(i);
%% Optimization Problem
disp('Starting optimization');

X_0_LMIs = diag(x_i);
% optimize(con,-log(det(X_max)))
options = sdpsettings('solver','sdpt3');
diagnostics =optimize(con,-log(det(X_0_LMIs)));

if diagnostics.problem ~= 0
    if diagnostics.problem ==4
        disp('maximum iteration')
    else
        flagerror = diagnostics.problem
        disp('something wrong in the optimization')
    end
end


Y_0_LMIs = value(Y_0_LMIs);
X_0_LMIs = value(X_0_LMIs);
Q_xij=value(Q_xij);
% X_max = value(X_max);

P = inv(X_0_LMIs);
K = Y_0_LMIs*P;
%eigP = eig(inv(X_0_LMIs));

%CHECK LMI
[A,B,C,D] = getABCD(constants.theta,constants);
A_cl = A+B*K;
C_cl = C+D*K;
M_check = A_cl'*P*A_cl+Q+C_cl'*blockDiagonalMatrix(Q,constants.pbar-1)*C_cl+K'*R*K-P;
E = eig(M_check)
if max(E)>1e-7
   error('Stability condition not satisfied for nominal') 
end
%%
%double check barP
%eig(value((A_cl*X)'*P*(A_cl*X)-P_xj{1,1}-P_xj{1,2}-P_xj{1,3}))
for j=1:3
test=max(eig(value((A_cl(j,:)*X)'*P(j,j)*(A_cl(j,:)*X)-P_xj{1,j})));
     if test>1e-7
       error(['error at:' num2str([num_state,test*1e7])]) 
     end
end
%%
for pred = 1:constants.pbar-1
    for num_state = 1:constants.nx
        C_temp=C_cl((pred-1)*constants.nx+num_state,:);
test=max(eig(value((C_temp*X)'*Q(j,j)*(C_temp*X)-Q_xij{pred,num_state})));
     if test>1e-7
       error(['error at:' num2str([pred,num_state,test*1e7])]) 
     end
    end
end 

end
