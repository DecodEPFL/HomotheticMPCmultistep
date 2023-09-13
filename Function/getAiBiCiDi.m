function [A,B,C,D] = getAiBiCiDi(theta,constants,pbar)
    if size(theta,2)>size(theta,1)
        theta = theta';
    end


    if pbar == constants.pbar
        C = 0;
        D = 0;
        A = theta(1:constants.nx)';
        B = flip(theta(constants.nx+1:end)');
    else
        A = 0;
        B = 0;
        C = theta(1:constants.nx)';
        D = [flip(theta(constants.nx+1:end)') zeros(1, constants.pbar-length(theta(constants.nx+1:end)'))];
    end
end