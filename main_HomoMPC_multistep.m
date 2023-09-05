clear all
close all
clc

load("offlineComputation.mat")


constants.n = size(multiStepModel.Abar,2);
constants.m = size(multiStepModel.Bbar,2);
constants.p = size(multiStepModel.Cbar,1);
constants.np = size(constants.theta,1);

constants.A = multiStepModel.Abar;
constants.B = multiStepModel.Bbar;
constants.C = multiStepModel.Cbar;
constants.D = multiStepModel.Dbar;

constants.Xf_homo = constants.X0;
constants.cU = constants.cU;
constants.F = constants.F;
constants.G = constants.G;
constants.Gbig = constants.Gbig;
constants.Fbig = constants.Fbig;
constants.X_0 = constants.X0;

constants.Hx = constants.Hx;
constants.hx = constants.hx;

constants.P = constants.P;
constants.K = constants.K;
constants.Q = constants.Q;
constants.R = constants.R;

constants.H_theta = constants.Htheta;
constants.h_theta = constants.htheta;

constants.N = 6;

% creation of the solver
MPC = create_multi_rate_mpc(constants);

% simulation time
Tsim = 15; 
% running simulation
x0 = [2;2;2];
[x_traj,u_traj, y_traj,lambda,alpha,flag]= simulation(MPC,x0,constants,Tsim);%simulate

XX = [];
for i = 1:Tsim
    XX =[XX;x_traj(:,i);y_traj(:,i)];
end
figure;
plot(XX(3:constants.nx:end));title('x3')

function yalmip_optimizer = create_multi_rate_mpc(constants)
% extract from constants the variable needed
v_0 = constants.v_0;
c_0 = constants.c_0;
q_t = constants.q_t;
Hx=constants.Hx;
hx=constants.hx;
F=constants.F;
Fbig=constants.Fbig;
G=constants.G;
N = constants.N;
A_0=constants.A_i(:,1:constants.n);
B_0=constants.B_i(:,1:constants.m);
C_0=constants.C_i(:,1:constants.n);
D_0=constants.D_i(:,1:constants.m);
wbarx = constants.wbarx;
wbary = constants.wbary;
pbar = constants.pbar;
H_theta = constants.H_theta;
h_theta = constants.h_theta;
nF = size(constants.F,1);
nx = constants.nx;

K=constants.K;
fbar = constants.fbar;
gbar = constants.gbar;
[A_bar,B_bar,C_bar,D_bar]=getABCD(constants.theta,constants);

% define symbolic optimization variables
x_t = sdpvar(constants.n,1,'full');%measured state
Z = sdpvar(constants.n,constants.N+1,'full');%nominal state z (center of tube)
V = sdpvar(constants.m,constants.N,'full'); 
Y = sdpvar(constants.p,constants.N,'full'); 
alpha = sdpvar(1,constants.N+1,'full');
%1. define Lambda as zero matrix
Indices = zeros(N,v_0,c_0+nF*(pbar-1),q_t);
%2. add sdpvar entries c_0
%first c_0 entries; which param: 
%JKIndices(:,:,1:c_0,1:(constants.pbar*constants.nu+constants.nx)*2)=ones(N,v_0,c_0,(constants.pbar*constants.nu+constants.nx)*2);
Indices(:,:,1:c_0,1:(constants.pbar*constants.nu+constants.nx)*2*nx)=ones(N,v_0,c_0,(constants.pbar*constants.nu+constants.nx)*2*constants.nx);
%JKLambda(:,:,c_0+1:end,(constants.pbar*constants.nu+constants.nx)*2+1:end)=0*Lambda(:,:,c_0+1:end,(constants.pbar*constants.nu+constants.nx)*2+1:end);
%3. add sdpvar entries for each F
%JKstart_index_Htheta=(constants.pbar*constants.nu+constants.nx)*2;

start_index_Htheta=(constants.pbar*constants.nu+constants.nx)*2*nx;
for p=constants.pbar-1:-1:1
    %index in H_theta corresponding to predictor p:
    index_HTheta=start_index_Htheta+1:start_index_Htheta+(p*constants.nu+constants.nx)*2*nx;
    start_index_Htheta=index_HTheta(end);
    %index in [X0;blkdiag(F)] corresponding to predictor p:
    index_F=c_0+nF*(p-1)+1:c_0+nF*p;
    Indices(:,:,index_F,index_HTheta)=ones(N,v_0,nF,length(index_HTheta));    
end

Lambda = sdpvar(N,v_0,c_0+nF*(pbar-1),q_t,'full');
% fix as 0 the variable in Lambda not used
Lambda=Lambda.*Indices;

%initialize cost
objective = 0;
constraints = [Hx*(x_t - Z(:,1)) - alpha(1)*hx<=0]; % initial state constraint 
for k = 1:N
    Y(:,k) = C_bar*Z(:,k)+D_bar*(K*Z(:,k)+V(:,k));
    % dynamics
    constraints = [constraints, F*Z(:,k) + alpha(k)*fbar<=ones(size(F,1),1)];
    constraints = [constraints, (G*K)*Z(:,k) + G*V(:,k) + alpha(k)*gbar<=ones(size(G,1),1)];
    for j=1:constants.v_0 %cycle through vertices of tube
        %compute vertices tube:
         x_kj=Z(:,k) + alpha(k)*constants.X_0.V(j,:)';
         u_kj= K*x_kj+V(:,k);    
        % computing e^kj and E^kj        
        E_kj= getBigE(x_kj,u_kj,constants);
        e_kj= [A_0*x_kj+B_0*u_kj-Z(:,k+1);
            C_0*x_kj+D_0*u_kj];
       
        phi =[alpha(k+1)*ones(c_0,1); ones(nF*(pbar-1),1)];
        Hpbar = blkdiag(Hx,Fbig);
        wbar = [wbarx;wbary];
        %(17d) Lambda_jk*h_theta+Hp*e_jk-phi_k<=-w_bar
        constraints = [constraints, reshape(Lambda(k,j,:,:),c_0+nF*(pbar-1),q_t)*h_theta+Hpbar*e_kj-phi+wbar<=0];
        %(17e): Hp*E_jk=Lambda_kj*H_theta
        constraints = [constraints, Hpbar*E_kj==reshape(Lambda(k,j,:,:),c_0+nF*(pbar-1),q_t)*H_theta];
        %Lambda>=0
            for l=1:c_0+nF*(pbar-1)
            constraints = [constraints, reshape(Lambda(k,j,l,:),1,q_t)>=0];
            end
    end
    %alpha>0
    objective = objective + 10*alpha(k);
    constraints = [constraints,alpha(k)>=0];
    objective = objective +stagecosts(Y(:,k),V(:,k),constants.Q,constants.R) ; %add stage cost
end

objective = objective + terminalcosts(Z(:,N+1),constants.P); %add terminal cost
objective = objective + sum(alpha(:))*1e-3;
% setup optimizer object
ops = sdpsettings('verbose',1,'solver','quadprog');%solver settings
yalmip_optimizer = optimizer(constraints,objective,ops,x_t,{objective,Z,V,Y,alpha,Lambda});%set solver
end


function [x_traj, u_traj, y_traj,lambda,alpha,flag] = simulation(yalmip_optimizer, x0, constants, Tsim)
%simulate an MPC given by "yalmip_optimizer"
    u_traj = [];
    x_traj = [x0];
    y_traj = [];
    %w_traj = [];
    u_sol_mpc = {}; 
    for k = 1:1:Tsim
        [res,flag] = yalmip_optimizer(x_traj(:,end)); % call mpc
        assert(flag==0||flag==3);%check for errors/infeasibility
        obj_sol_mpc{k} = res{1}; %optimal cost
        x_sol_mpc{k} = res{2};%get nominal state
        u_sol_mpc{k} = res{3};%get nominal input
        lambda{k} = res{6};%get lambda
        alpha{k} = res{5};%get alpha
        
        %w_traj(:,end+1) = constants.W.V(1,:)';%simply take one vertex for now
        u_traj(:,end+1) = u_sol_mpc{k}(:,1)+constants.K*x_traj(:,k);%apply input
        x_traj(:,end+1) = constants.A * x_traj(:,end) +constants.B * u_traj(:,end);%simulate one step
        y_traj(:,end+1) = constants.C * x_traj(:,end) +constants.D * u_traj(:,end);
    end
    x_traj(:,end)=[];%remove last state (to have trajectories of length Tsim)
end 

function [A,B,C,D] = getABCD(theta,constants)
    %A,B for given parameter theta
    temp=kron([1;theta],eye(constants.n));
    A=constants.A_i*temp;
    C=constants.C_i*temp;
    temp=kron([1;theta],eye(constants.m));
    B=constants.B_i*temp;
    D=constants.D_i*temp;
end
function [out] = getBigE(x,u,constants) 
    temp_x=kron(eye(constants.np),x);
    temp_u=kron(eye(constants.np),u);
    x = constants.A_i(:,constants.n+1:end)*temp_x+constants.B_i(:,constants.m+1:end)*temp_u;%leave out A_0,B_0,C_0,D_0
    y = constants.C_i(:,constants.n+1:end)*temp_x+constants.D_i(:,constants.m+1:end)*temp_u;%leave out A_0,B_0,C_0,D_0
    out =[x;y];
end


function cost = stagecosts(y,u,Q,R)  
    cost = y'*Q*y + u'*R*u;
end

function cost = terminalcosts(x,P)
   cost = x'*P*x;
end