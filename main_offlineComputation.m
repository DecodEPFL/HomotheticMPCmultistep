bclear all
close all
clc

addpath('Function')  
load IdentifiedModelandData.mat

ppred = Bound.ppred;
pbar = ppred;
nx = RealSys.nx;

nw = 1;%RealSys.nw;
M = RealSys.M;
input = data.u;
state = data.xx;

constants.nx = nx;
constants.nu = RealSys.nu;
constants.pbar = pbar;

% test of the obtained multistep predictors
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

% input and distrubance
j = 1;
for i = 0:ppred:length(input)-ppred
    Ubig{j} = input(i+1:i+ppred);
    Wbig{j} = data.w(i+1:i+ppred);
    j = j+1;
end

% building the true matrices
for i = 1:nx
    Abar(i,:) = theta{i,ppred}(1:nx);
    Bbar(i,:) = theta{i,ppred}(nx+1:end);
    AbarTRUE(i,:) = thetaTRUE{i,ppred}(1:nx);
    BbarTRUE(i,:) = thetaTRUE{i,ppred}(nx+1:end);
end

% building the big matrices with the identified parameters y = phi^T*theta
for j = 1:ppred-1
    for i = 1:nx
        Cbar((j-1)*nx+i,:) = theta{i,j}(1:nx);
        Dbar((j-1)*nx+i,1:length(theta{i,j}(nx+1:end))) = theta{i,j}(nx+1:end);
        CbarTRUE((j-1)*nx+i,:) = thetaTRUE{i,j}(1:nx);
        DbarTRUE((j-1)*nx+i,1:length(thetaTRUE{i,j}(nx+1:end))) = thetaTRUE{i,j}(nx+1:end);
    end
    iter = 0;
    for index = j:-1:1
        NbarTRUE((j-1)*nx+1:(j-1)*nx+nx,index) = RealSys.A^(iter)*RealSys.M;
        iter = iter+1;
    end
end

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

for jp = 1:Nbig
    X1big = Abar*X0big+Bbar*Ubig{jp}'+Mbar*Wbig{jp}';
    Ybig = Cbar*X0big + Dbar*Ubig{jp}' + Nbar*Wbig{jp}';
    X0big = X1big;
    YY = [YY; Ybig; X1big];
end

% Plot of the data
figure;plot(state(:,1));hold on;plot(YY(1:nx:end-nx));  legend('data','estimated');title('x1')
figure;plot(state(:,2));hold on;plot(YY(2:nx:end-nx));  legend('data','estimated');title('x2')
figure;plot(state(:,3));hold on;plot(YY(3:nx:end-nx));  legend('data','estimated');title('x3')

% Test of the obtained multi-step model
X0big = state(1,:)';
YY = [];
YY = [YY; X0big];
for jp = 1:Nbig
    X1big = AbarTRUE*X0big+BbarTRUE*Ubig{jp}'+MbarTRUE*Wbig{jp}';
    Ybig = CbarTRUE*X0big + DbarTRUE*Ubig{jp}' + NbarTRUE*Wbig{jp}';
    X0big = X1big;
    YY = [YY; Ybig; X1big];
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

%theta =[x1pbar, x2pbar, x3pbar, ...,...x11,x21,x31]

num_row = 0; num_rowP = 0;
num_col = 0; num_colP = 0;
for j = 1:ppred
    for i = 1:nx
        if FPS{i,j}.isEmptySet
            min(FPS{i,j}.A*thetaTRUE{i,j}'<=FPS{i,j}.b)
            FPS{i,j}
            FPS{i,j}.isEmptySet
        end
        % outer approximation of the FPS
        Ftemp = outerApproxBox(FPS{i,j}.A,FPS{i,j}.b);
        
        FPS{i,j}= Ftemp;
        num_row = num_row + size(FPS{i,j}.A,1);
        num_col = num_col + size(FPS{i,j}.A,2);
        if j == ppred
            num_rowP = num_rowP + size(FPS{i,j}.A,1);
            num_colP = num_colP + size(FPS{i,j}.A,2);
        end
    end
end

% building overall FPS and stacking parameters
FPSA = zeros(num_row,num_col);
FPSAP = zeros(num_rowP,num_colP);
FPSb = []; FPSbP = [];
thetaPbar = [];
for j = 1:constants.pbar
    for i = constants.nx:-1:1
        if FPS{i,j}.A*RealSys.thetaTRUE{j,i}'>FPS{i,j}.b
            stop =1;
        end
        if FPS{i,j}.A*theta{i,j}'>FPS{i,j}.b
            stop =1;
        end
        thetaTRUEP = [RealSys.thetaTRUE{j,i}';thetaTRUEP];
        thetaP = [theta{i,j}';thetaP];
        matrA = FPS{i,j}.A;
        colA = size(matrA,2); rowA = size(matrA,1);
        num_row = num_row-rowA; num_col = num_col-colA;
        FPSA(num_row+1:num_row+rowA,num_col+1:num_col+colA)=matrA;
        FPSb=[FPS{i,j}.b;FPSb];
        if j ==pbar
            thetaPbar = [theta{i,j}';thetaPbar];
            matrA = FPS{i,j}.A;
            colA = size(matrA,2); rowA = size(matrA,1);
            num_rowP = num_rowP-rowA; num_colP = num_colP-colA;
            FPSAP(num_rowP+1:num_rowP+rowA,num_colP+1:num_colP+colA)=matrA;
            FPSbP=[FPS{i,j}.b;FPSbP];
        end
    end
end

if max(FPSA*thetaP-FPSb)<1e-8
    disp('theta p not belongs to FPS')
end
if max(FPSA*thetaTRUEP-FPSb)<1e-8
    disp('theta p TRUE not belongs to FPS')
end

FPS_big = Polyhedron('A',FPSA,'b',FPSb); 
FPSP = Polyhedron('A',FPSAP,'b',FPSbP); 

constants.thetapbar = thetaPbar;
constants.theta = thetaP;
constants.thetaTRUE = thetaTRUEP;
constants.np = size(thetaP,1);

nx = constants.nx;

%% build Ai Bi Ci Di matrices
[A_i,B_i,C_i,D_i,A_ip,B_ip] = buildAiBiCiDi(constants,multiStepModel,thetaP,thetaPbar,theta);
constants.A_i=A_i;
constants.B_i=B_i;
constants.C_i=C_i;
constants.D_i=D_i;
constants.A_ip=A_ip;
constants.B_ip=B_ip;
nu_ms = constants.nu*constants.pbar;
constants.nu_ms = nu_ms;

% DEFINITION OF THE SETS
%constraints on input
umax=constants.umax;
umin = constants.umin;

Au=[eye(pbar)./umax;
    -eye(pbar)./(-umin)];
bu=[ones(pbar,1);
    ones(pbar,1)];

cU=Polyhedron('A',Au,'b',bu);
constants.cU = cU;

constants.G = Au;


% constraints on state
xmax = constants.xmax;
xmin = constants.xmin;
Ax=[eye(nx); 
    -eye(nx)];
bx=[ones(nx,1)*xmax;
    ones(nx,1)*(-xmin)];
cX=Polyhedron('A',Ax,'b',bx);


constants.F = Ax./bx;



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
%% Design of low complexity X0
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
X0=Polyhedron([inv(V);-inv(V)],ones(2*nx,1));


constants.Hx = X0.A;
constants.hx = X0.b;
constants.X0 = X0;
constants.c_0 = length(X0.b(:,1));
constants.v_0 = length(X0.V(:,1));
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

Hw = [eye(pbar);-eye(pbar)];
hw = [data.wmax*ones(2*pbar,1)];
W = Polyhedron(Hw,hw);
Wx = MbarTRUE*W;
constants.Wx = Wx;

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
Wy = NbarTRUE*W;
constants.Wy = Wy;

% definition of Fbig
Fbig = blockDiagonalMatrix(constants.F, pbar-1);
constants.Fbig = Fbig;
% definition of Gbig
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
constants.alphamax = 0.034611662920428;


save('offlineComputation4','RealSys','constants','data','multiStepModel','Bound')


% Function to outer-approximate a polytopic set with a box
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
        box = Polyhedron('A',[eye(p);-eye(p)],'b',[upper';-lower']);
    else

        box = box;
    end
end

