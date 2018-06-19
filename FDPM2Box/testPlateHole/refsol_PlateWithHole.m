function [sol] = refsol_PlateWithHole(E, nu, a, Tx)

mu = E / (2*(1+nu));
k = 3 - 4*nu; % plane strain
% k = (3-nu)/(1+nu); % plane strain

a2 = a^2;
a4 = a^4;

function [ur,ut] = sol_disp_rtheta(r, theta)
    c2t = cos(theta.*2);
    s2t = sin(theta.*2);
    
    rinv = 1 ./ r;
    r3 = r.^3;
    
    ur = r.*((k-1)/2+c2t) + a2.*rinv.*(1+(1+k).*c2t) - a4./r3.*c2t;
    ur = Tx/(4*mu) .* ur;
    
    ut = ((1-k)*a2.*rinv - r - a4./r3) .* s2t;
    ut = Tx/(4*mu) .* ut;
end

function [ux,uy] = sol_disp_xy(x, y)
    r = hypot(x, y);
    theta = atan2(y, x);
    
    [ur,ut] = sol_disp_rtheta(r, theta);
    
    st = sin(theta);
    ct = cos(theta);
    ux = ur.*ct - ut.*st;
    uy = ur.*st + ut.*ct;
end

function [ux] = sol_ux(x, y)
    [ux,~] = sol_disp_xy(x,y);
end
function [uy] = sol_uy(x, y)
    [~,uy] = sol_disp_xy(x,y);
end

function [sxx,syy,sxy] = sol_stress(x, y)
    r = hypot(x, y);
    theta = atan2(y, x);
    
    r2 = r.^2;
    r4 = r.^4;
    
    c2t = cos(theta.*2);
    s2t = sin(theta.*2);
    c4t = cos(theta.*4);
    s4t = sin(theta.*4);
    
    sxx = Tx .* (1 - a2./r2.*(1.5*c2t+c4t) + 1.5*a4./r4.*c4t);
    syy = -Tx .* (a2./r2.*(0.5*c2t-c4t) + 1.5*a4./r4.*c4t);
    sxy = -Tx .* (a2./r2.*(0.5*s2t+s4t) - 1.5*a4./r4.*s4t);
end

function [sxx] = sol_sxx(x, y)
    [sxx,~,~] = sol_stress(x, y);
end
function [syy] = sol_syy(x, y)
    [~,syy,~] = sol_stress(x, y);
end
function [sxy] = sol_sxy(x, y)
    [~,~,sxy] = sol_stress(x, y);
end

% return a closure
sol = struct();
sol.disp_rt = @sol_disp_rtheta;
sol.disp_xy = @sol_disp_xy;
sol.stress = @sol_stress;
sol.ux = @sol_ux;
sol.uy = @sol_uy;
sol.sxx = @sol_sxx;
sol.syy = @sol_syy;
sol.sxy = @sol_sxy;





return
end




