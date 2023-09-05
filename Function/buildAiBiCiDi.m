function [A_i,B_i,C_i,D_i,A_ip,B_ip] = buildAiBiCiDi(constants,multiStepModel,thetaP,thetaPbar,theta)
nx = constants.nx;

A_i =[];
A_0 = zeros(constants.nx,constants.nx);
A_i = [A_i A_0];

for p = constants.pbar:-1:1
    for state = 1:nx
            for j = 1:nx
                var = state;
%                 if state == 1
%                     var = 3;
%                 elseif state == 3
%                     var = 1;
%                 else 
%                     var =2;
%                 end
                if p <= constants.pbar-1
                    A_temp = zeros(nx,nx);
                    A_i=[A_i A_temp];
                else
                    A_temp = zeros(nx,nx);
                    A_temp(var,j)=1;
                    A_i=[A_i A_temp];
                end
            end
            for m = 1:p
                A_i=[A_i zeros(nx,nx)];
            end
    end
end

constants.A_i=A_i;
temp=kron([1;thetaP],eye(constants.nx));
Abar = A_i*temp
AbarTRUE = multiStepModel.Abar



A_ip =[];
A_ip = [A_ip A_0];
A_ipi ={};

p = constants.pbar
for state = 1:nx
        for j = 1:nx
            var = state;
%             if state == 1
%                 var = 3;
%             elseif state == 3
%                 var = 1;
%             else 
%                 var =2;
%             end
            if p <= constants.pbar-1
                A_temp = zeros(nx,nx);
                A_ip=[A_ip A_temp];
            else
                A_temp = zeros(nx,nx);
                A_temp(var,j)=1;
                A_ip=[A_ip A_temp];
            end
        end
        for m = 1:p
            A_ip=[A_ip zeros(nx,nx)];
        end
end

constants.A_ip=A_ip;
temp=kron([1;thetaPbar],eye(constants.nx));
Abar = A_ip*temp
AbarTRUE = multiStepModel.Abar

for i = 1:constants.nx
    p = theta{i,constants.pbar};
    [Api,Bpi,Cpi,Dpi ] = getAiBiCiDi(p,constants,constants.pbar);
    Api'
end

% B_i
nu_ms = constants.nu*constants.pbar;
constants.nu_ms = nu_ms;
nx = constants.nx;


B_0 = zeros(nx,nu_ms);
B_i = [];
B_i = [B_i B_0];
for p = constants.pbar:-1:1
    for state = 1:nx
            for j = 1:nx
                B_i=[B_i zeros(nx,nu_ms)];
            end
            if p == constants.pbar
                for m = 1:p
                    B_temp = zeros(nx,nu_ms);
                          var = state;
%                         if state == 1
%                             var =3;
%                         elseif state == 3
%                             var =1;
%                         else 
%                             var = 2;
%                         end
                    B_temp(var,m)=1;
                    B_i=[B_i B_temp];
                end
            else
                for m = 1:p
                B_temp =zeros(nx,nu_ms);
                B_i=[B_i B_temp];
                end
            end
    end
end

temp=kron([1;thetaP],eye(constants.nu*constants.pbar));
Bbar = B_i*temp
BbarTRUE=multiStepModel.Bbar
constants.B_i=B_i;

B_0 = zeros(nx,nu_ms);
B_ip = [];
B_ip = [B_ip B_0];
p = constants.pbar
for state = 1:nx
        for j = 1:nx
            B_ip=[B_ip zeros(nx,nu_ms)];
        end
        if p == constants.pbar
            for m = 1:p
                B_temp = zeros(nx,nu_ms);
                      var = state;
%                     if state == 1
%                         var =3;
%                     elseif state == 3
%                         var =1;
%                     else 
%                         var = 2;
%                     end
                B_temp(var,m)=1;
                B_ip=[B_ip B_temp];
            end
        else
            for m = 1:p
            B_temp =zeros(nx,nu_ms);
            B_ip=[B_ip B_temp];
            end
        end
end

temp=kron([1;thetaPbar],eye(constants.nu*constants.pbar));
Bbar = B_ip*temp
BbarTRUE=multiStepModel.Bbar
constants.B_ip=B_ip;
for i = 1:constants.nx
    p = theta{i,constants.pbar};
    [Api,Bpi,Cpi,Dpi ] = getAiBiCiDi(p,constants,constants.pbar);
    Bpi'
end


% C_i
ny_ms = (constants.pbar-1)*nx;
C_0 = zeros(ny_ms,nx);
C_i = [];
C_i = [C_i C_0];
for p = constants.pbar:-1:1
    for state = 1:nx
            for j = 1:nx
                  var = state;
%                 if state == 1
%                     var = 3;
%                 elseif state == 3
%                     var = 1;
%                 else 
%                     var =2;
%                 end
                if p <= constants.pbar-1
                C_temp = zeros(ny_ms,nx);
                C_temp((p-1)*nx+var,j)=1;
                C_i=[C_i C_temp];
                else
                    C_temp = zeros(ny_ms,nx);
                    C_i=[C_i C_temp];
                end
            end
            for m = 1:p
                C_i=[C_i zeros(ny_ms,nx)];
            end
    end
end
constants.C_i=C_i;

temp=kron([1;thetaP],eye(constants.nx));
Cbar = C_i*temp
CbarTRUE = multiStepModel.Cbar


CpiT = [];
for pred = 1:constants.pbar-1
    for i = 1:constants.nx
         p = theta{i,pred};
        [Api,Bpi,Cpi,Dpi ] = getAiBiCiDi(p,constants,pred);
        CpiT = [CpiT;Cpi];
    end
end
CpiT

% D_i
D_0 = zeros(ny_ms,nu_ms);
D_i = [];
D_i = [D_i D_0];
for p = constants.pbar:-1:1
    for state = 1:nx
            for j = 1:nx
                D_i=[D_i zeros(ny_ms,nu_ms)];
            end
            for m = 1:p
                if p <= constants.pbar-1
                D_temp = zeros(ny_ms,nu_ms);
                  var = state;
%                 if state == 1
%                     var =3;
%                 elseif state == 3
%                     var =1;
%                 else 
%                     var = 2;
%                 end
                    D_temp((p-1)*nx+var,m)=1;
                    D_i=[D_i D_temp];
                else
                    D_temp = zeros(ny_ms,nu_ms);
                    D_i=[D_i D_temp];
                end
            end
    end
end
constants.D_i = D_i;
temp=kron([1;thetaP],eye(nu_ms));
Dbar = D_i*temp
DbarTRUE=multiStepModel.Dbar

DpiT = [];
for pred = 1:constants.pbar-1
    for i = 1:constants.nx
         p = theta{i,pred};
        [Api,Bpi,Cpi,Dpi ] = getAiBiCiDi(p,constants,pred);
        DpiT = [DpiT;Dpi];
    end
end
DpiT
end