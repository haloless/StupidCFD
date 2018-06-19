function [xs,ys, dxdphi,dydphi] = ParamCurve(phi,phi0,C,M)
%

sphi = sin(phi);
cphi = cos(phi);

tol = 1.0e-6;

if abs(M) > tol
    % integral solution
    
    coef = 1 / (-2*M);
    % coef = 1 / (2*M);
    
    
    % sign of solution
    % concave bridge, y is minus, x is plus
    ssy = -1;
    ssx = +1;
    
    % ssy = +1;
    % ssx = +1;
    
    ss = zeros(size(phi));
    ss(:) = ssy;
    % ss(1) = +1;
    % ss(end) = +1;
    ss(:) = +1;
    
    dd = sqrt(cphi.^2 + 4*C*M);
    %
    ys = coef * (cphi + ss.*dd);
    
    %
    k2 = 1 / (1+4*C*M);
    k = sqrt(k2);
    
    %
    ee = ellipticE(phi,k2) - ellipticE(phi0,k2);
    ff = ellipticF(phi,k2) - ellipticF(phi0,k2);
    
    xs = coef * (-sphi+sin(phi0) - ss.*(1/k*ee - 4*C*M*k*ff));
    
    % calculate derivative dx/dphi
    if nargout > 2
        dxdphi = coef * (-cphi - ss .* cphi.^2 ./ sqrt(cphi.^2+4*C*M));
        dydphi = 0;
    end
    
else
    % M=0, degenerate solution
    disp('degenerate');
    ys = C ./ cphi;
    
    xs = C * 2 * (atanh(tan(phi/2)) - atanh(tan(phi0/2)));
    
    
    % calculate derivative dx/dphi
    if nargout > 2
        dxdphi = C ./ cphi;
        dydphi = 0;
    end
end




return
end


