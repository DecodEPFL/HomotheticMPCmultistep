% Clear workspace, close figures, and clear command window
clear all
close all
clc

% Add path to function folder
addpath('Function')  

% Load data from a .mat file
load IdentifiedModelandData.mat

%% constraints
% Define input constraints
umax=10;
umin = -10;
constants.umax = umax;
constants.umin = umin;

% Define state constraints
xmax = 10;
xmin = [-10; -10; -10];
constants.xmax = xmax;
constants.xmin = xmin;

%% Get necessary parameters and matrices from loaded data
ppred = Bound.ppred;
pbar = ppred;
nx = RealSys.nx;

nw =RealSys.nw;
M = RealSys.M;
input = data.u;
state = data.xx;

constants.nx = nx;
constants.nu = RealSys.nu;
constants.pbar = pbar;

Bound.tau = Bound.tau_cons;

%% multistep predictors
thetaTRUE = RealSys.thetaTRUE';
AbarTRUE = zeros(nx,nx);
BbarTRUE = zeros(nx,ppred);
CbarTRUE = zeros(nx*(ppred-1),nx);
DbarTRUE = zeros(nx*(ppred-1),ppred);
MbarTRUE = zeros(nx,ppred);
NbarTRUE = zeros(nx*(ppred-1),ppred);

Abar = zeros(nx,nx);
Bbar = zeros(nx,ppred);
Mbar = zeros(nx,ppred);
Mbar(:,ppred) = M;

Cbar = zeros(nx*(ppred-1),nx);
Dbar = zeros(nx*(ppred-1),ppred);
Nbar = blockDiagonal(M,ppred-1);

Nbar = [Nbar zeros(nx*(ppred-1),1)];

j = 1;
% definition of exogenous signals
for i = 0:ppred:length(input)-ppred
    Ubig{j} = input(i+1:i+ppred);
    Wbig{j} = data.w(i+1:i+ppred);
    j = j+1;
end

% Abar Bbar
for i = 1:nx
    Abar(i,:) = theta{i,ppred}(1:nx);
    Bbar(i,:) = flip(theta{i,ppred}(nx+1:end));
    AbarTRUE(i,:) = thetaTRUE{i,ppred}(1:nx);
    BbarTRUE(i,:) = flip(thetaTRUE{i,ppred}(nx+1:end));
end
% Dbar Cbar Nbar
for j = 1:ppred-1
    for i = 1:nx
        Cbar((j-1)*nx+i,:) = theta{i,j}(1:nx);
        Dbar((j-1)*nx+i,1:length(theta{i,j}(nx+1:end))) = flip(theta{i,j}(nx+1:end));
        CbarTRUE((j-1)*nx+i,:) = thetaTRUE{i,j}(1:nx);
        DbarTRUE((j-1)*nx+i,1:length(thetaTRUE{i,j}(nx+1:end))) = flip(thetaTRUE{i,j}(nx+1:end));
    end
    iter = 0;
    for index = j:-1:1
        NbarTRUE((j-1)*nx+1:(j-1)*nx+nx,index) = RealSys.A^(iter)*RealSys.M;
        iter = iter+1;
    end
end

% definition Mbar
col = 1;
for j = ppred-1:-1:1
    MbarTRUE(:,col) = RealSys.A^(j)*RealSys.M;
    col = col +1;
end
MbarTRUE(:,ppred) =RealSys.M;

X0big = state(1,:)';
Nbig = length(Ubig);
YY = [];
YY = [YY; X0big];

% Test identified multistep
for jp = 1:Nbig
    X1big = Abar*X0big+Bbar*Ubig{jp}';%MbarTRUE*Wbig{jp}';
    Ybig = Cbar*X0big + Dbar*Ubig{jp}';% NbarTRUE*Wbig{jp}';
    X0big = X1big;
    YY = [YY; Ybig; X1big];
end


figure;plot(state(:,1));hold on;plot(YY(1:nx:end-nx));  legend('data','estimated');title('x1')
figure;plot(state(:,2));hold on;plot(YY(2:nx:end-nx));  legend('data','estimated');title('x2')
figure;plot(state(:,3));hold on;plot(YY(3:nx:end-nx));  legend('data','estimated');title('x3')


X0big = state(1,:)';
YY = [];
YY = [YY; X0big];
for jp = 1:Nbig
    X1big = AbarTRUE*X0big + BbarTRUE*Ubig{jp}' + MbarTRUE*Wbig{jp}';
    Ybig = CbarTRUE*X0big + DbarTRUE*Ubig{jp}' + NbarTRUE*Wbig{jp}';
    X0big = X1big;
    YY = [YY; Ybig];
    YY = [YY; X1big];
end

figure;plot(state(:,1));hold on;plot(YY(1:nx:end-nx));title('trueparam multistep');legend('data','estimated')

XX=[];
x0 = zeros(3,1);
for i = 1:size(input,2)
    x1 = RealSys.A*x0+RealSys.B*input(i)+RealSys.M*data.w(i);
    XX = [XX;x0];
    x0 = x1;
end
figure;plot(state(:,1));hold on;plot(XX(1:nx:end-nx));title('trueparam statespace');legend('data','estimated')


%% SAVE THE MULTISTEP MODEL

multiStepModel.Abar = Abar;
multiStepModel.Bbar = Bbar;
multiStepModel.Cbar = Cbar;
multiStepModel.Dbar = Dbar;
multiStepModel.Mbar = Mbar;
multiStepModel.Nbar = Nbar;

multiStepModel.AbarTRUE = AbarTRUE;
multiStepModel.BbarTRUE = BbarTRUE;
multiStepModel.CbarTRUE = CbarTRUE;
multiStepModel.DbarTRUE = DbarTRUE;
multiStepModel.MbarTRUE = MbarTRUE;
multiStepModel.NbarTRUE = NbarTRUE;
multiStepModel.theta = theta;

%% Compute A_i B_i C_i matrices and FPS

thetaTRUEP = [];
thetaP = [];

FPSA =[];
FPSb =[];
for j = 1:constants.pbar
    for i = constants.nx:-1:1
        if FPS{i,j}.isEmptySet
            min(FPS{i,j}.A*thetaTRUE{i,j}'<=FPS{i,j}.b)
            error('FPS empty')
        end
        if max(FPS{i,j}.A*RealSys.thetaTRUE{j,i}'-FPS{i,j}.b)>1e-7
            error('FPS does not contain true param')
        end
        FPS{i,j} = outerApproxBox(FPS{i,j}.A,FPS{i,j}.b);
        thetaTRUEP = [RealSys.thetaTRUE{j,i}';thetaTRUEP];
        thetaP = [theta{i,j}';thetaP];
        matrA = FPS{i,j}.A;
        FPSA = blkdiag(matrA,FPSA);
        FPSb=[FPS{i,j}.b;FPSb];
        if max(FPSA*thetaP-FPSb)>1e-7
            error('parameters not in FPS')
        end
    end
end

if max(FPSA*thetaP-FPSb)>1e-8
    disp('theta p not belongs to FPS')
end
FPS_big = Polyhedron('A',FPSA,'b',FPSb); 
if FPS_big.isEmptySet
    error('FPS is empty')
end
constants.theta = thetaP;
constants.thetaTRUE = thetaTRUEP;
constants.np = size(thetaP,1);

nx = constants.nx;

%build AiBiCiDi
[A_i,B_i,C_i,D_i] = buildAiBiCiDi(constants,multiStepModel,thetaP,theta);
constants.A_i=A_i;
constants.B_i=B_i;
constants.C_i=C_i;
constants.D_i=D_i;
nu_ms = constants.nu*constants.pbar;
constants.nu_ms = nu_ms;

%% DEFINITION OF THE SETS
%constraints on input
umax = constants.umax;
umin = constants.umin;

constants.G=[eye(constants.pbar)./umax;
            -eye(constants.pbar)./(-umin)];
bu=[ones(constants.pbar,1);
    ones(constants.pbar,1)];

constants.cU = Polyhedron(constants.G,bu);


% constraints on state
xmax = constants.xmax;
xmin = constants.xmin;
constants.F=[eye(nx)./xmax; 
    -eye(nx)./(-xmin)];
bx=[ones(nx,1);
    ones(nx,1)];
cX=Polyhedron(constants.F,bx);


% FPS
constants.FPS = FPS;
constants.Theta = FPS_big;
constants.Htheta = FPS_big.A;
constants.htheta = FPS_big.b;

%%
%%%% AUXILIARY CONTROL LAW & TERMINAL COST

ps_k = 10;    
pu_k = 1;

Q = blkdiag(ps_k*eye(nx));
QY = blockDiagonal(Q,ppred-1);
R = pu_k*eye(pbar);

constants.Q = QY;
constants.R = R;

[P,K] = getPK_LMIs_D(Q,R,constants);

constants.P = P;
constants.K = K;
%%
Hw = [eye(pbar);-eye(pbar)];
hw = [data.wmax*ones(2*pbar,1)];
W = Polyhedron(Hw,hw);
Wx = MbarTRUE*W; %or use estimate lambda
constants.Wx = Wx;
Wy = NbarTRUE*W;%or use estimate lambda
constants.Wy = Wy;

A_cl = Abar+Bbar*K;
[V,lambda]=eig(A_cl);

if max(abs(eig(A_cl))>1)
    disp('Unstable system')
end

if min(min(imag(V))==[0 0 0])
    V = V;
elseif min(imag(V(:,1))==0)
    V=([V(:,1),real(V(:,2)),imag(V(:,2))]);
elseif min(imag(V(:,3))==0)
    V=([V(:,3),real(V(:,2)),imag(V(:,2))]);
end

X0 = Polyhedron([inv(V);-inv(V)],ones(2*nx,1));

constants.Hx = [inv(V);-inv(V)];
constants.hx = ones(2*nx,1);
constants.X0 = X0;

constants.c_0 = length(constants.X0.b(:,1));
constants.v_0 = length(constants.X0.V(:,1));
constants.q_t = length(FPS_big.b(:,1));  


% Computation of $\bar{f}$
fbar = [];
for i=1:length(constants.F(:,1))
    [~,fval] = linprog(-constants.F(i,:), constants.Hx, constants.hx,[],[],[],[],[],[]);
    if isempty(fval)
        disp('error in computing fbar')
    end
    fbar = [fbar; -fval];
end

% Computation of $\bar{g}$
gbar = [];
GK = (constants.G*constants.K);
for i=1:length(GK(:,1))
    [~,fval] = linprog(-GK(i,:), constants.Hx, constants.hx,[],[],[],[],[],[]);
    if isempty(fval)
        disp('error in computing gbar')
    end
    gbar = [gbar; -fval];
end

constants.fbar = fbar;
constants.gbar = gbar;


% Computation of $\bar{w}x$
wbarx = [];
for i=1: length(constants.Hx(:,1))
    [~,fval] = linprog(-constants.Hx(i,:), Wx.A, Wx.b,[],[],[],[],[],[]);
    if isempty(fval)
        disp('error in computing wbarx');
    end
    wbarx = [wbarx; -fval];
end

constants.wbarx = wbarx;

Fbig = blockDiagonalMatrix(constants.F, pbar-1);
constants.Fbig = Fbig;

Gbig = blockDiagonalMatrix(constants.G, pbar-1);
constants.Gbig = Gbig;

% Computation of $\bar{w}y$
wbary = [];
for i=1: length(constants.Fbig(:,1))
    [~,fval] = linprog(-constants.Fbig(i,:), Wy.A, Wy.b,Wy.Ae,Wy.be,[],[],[],[]);
    if isempty(fval)
        disp('error in computing wbary')
    end
    wbary = [wbary; -fval];
end

constants.wbary = wbary;
constants.alphamax = 3.34;


save('offlineComputation','RealSys','constants','data','multiStepModel','Bound')




function box = outerApproxBox(A,b)
    p = size(A,2);
    options=optimset('Display','none','tolcon',1e-8);
    for dim = 1 : p
        [~,lower(dim)]=linprog([zeros(dim-1,1);1;zeros(p-dim,1)],A,b,[],[],[],[],options);
         [~,upper(dim)]=linprog(-[zeros(dim-1,1);1;zeros(p-dim,1)],A,b,[],[],[],[],options);
        upper(dim) = -upper(dim);
    end

    A=[eye(p);-eye(p)];
    b=[upper';-lower'];
    boxT = polytope(A,b);
    box = Polyhedron('A',A,'b',b,'V',extreme(boxT));
    if box.isEmptySet
        %box = Polyhedron('Ae',eye(p),'be',[lower]);
        box = Polyhedron('A',[eye(p);-eye(p)],'b',[upper';-lower']);
    else

        box = box;
    end
end

