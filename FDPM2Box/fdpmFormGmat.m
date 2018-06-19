function [G] = fdpmFormGmat(np, dNx,dNy, prob_type,x,y)
%fdpmFormGmat
% G is 4-component full strain-displacement matrix
% [11,21,12,22]
% [Nx,0; 0,Nx; Ny,0; 0,Ny]

ndim = 2;

if prob_type == 1
    % plane-strain
    G = zeros(4,np*ndim);
elseif prob_type == 2
    % axisymmetric
    G = zeros(5,np*ndim);
end

% 11: [Nx,0]
G(1,1:ndim:end) = dNx;

% 21: [0,Nx]
G(2,2:ndim:end) = dNx;

% 12: [Ny,0]
G(3,1:ndim:end) = dNy;

% 22: [0,Ny]
G(4,2:ndim:end) = dNy;

% RZ coord
if prob_type == 2
    rr = x(end);
    if rr == 0
        % G(5,end-1) = 0.0;
        G(5,1:ndim:end) = dNx;
    else
        G(5,end-1) = 1.0 / rr;
    end
end

return
end

