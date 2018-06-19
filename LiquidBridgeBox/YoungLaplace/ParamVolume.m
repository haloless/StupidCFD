function [vol] = ParamVolume(phi1,phi2, C,M)
    
    fun = @(phi) volfunc(phi, phi1, C,M);
    
    if 1
        vol = integral(fun, phi1,phi2, 'AbsTol',1.0e-10,'RelTol',0.0);
    else
        % ndiv = 21;
        % ndiv = 51;
        % ndiv = 101;
        ndiv = 201;
        % ndiv = 501;
        phis = linspace(phi1,phi2,ndiv);
        % dvol = volfunc(phis, phi1, C,M);
        % vol = trapz(phis,dvol);
        [x,y] = ParamCurve(phis, phi1, C,M);
        vol = trapz(x,pi*y.^2);
    end
    
return
end


function [dvol] = volfunc(phi, phi0, C,M)
    
    [x,y,dxdphi,~] = ParamCurve(phi, phi0, C,M);
    
    dvol = pi * y.^2 .* dxdphi;
    % dvol = pi * y.^2;
    
return
end


