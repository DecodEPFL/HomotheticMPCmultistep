function [A,B,C,D] = getABCD(theta,constants)
    %A,B for given parameter theta
    temp=kron([1;theta],eye(constants.nx));
    A=constants.A_i*temp;
    C=constants.C_i*temp;
    temp=kron([1;theta],eye(constants.nu_ms));
    B=constants.B_i*temp;
    D=constants.D_i*temp;
end