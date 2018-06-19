function [B] = fdpmFormBmat(np,dNx,dNy, prob_type,x,y)
%fdpmFormBmat
% B is the 3-component standard strain-displacement
% [Nx,0; 0,Ny; Ny,Nx]

ndim = 2;

if ~exist('prob_type', 'var')
    prob_type = 1;
end

if prob_type == 1
    % plane-strain
    B = zeros(3,np*ndim);
elseif prob_type == 2
    % axisymmetric
    B = zeros(4,np*ndim);
end

% 11: [Nx, 0]
B(1,1:ndim:end) = dNx;

% 22: [0, Ny]
B(2,2:ndim:end) = dNy;

% 12: [Ny,Nx]
B(3,1:ndim:end) = dNy;
B(3,2:ndim:end) = dNx;

% RZ coord
if prob_type == 2
    % in this case, the hoop deformation is u/r
    % NOTE the center particle is always put in the last point
    rr = x(end);
    if rr == 0
        % B(4,end-1) = 0;
        B(4,1:ndim:end) = dNx;
    else
        B(4,end-1) = 1.0 / rr;
    end
end


return
end




