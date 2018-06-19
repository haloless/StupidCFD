function [sol] = refsol_PresLoadHalfPlane(E,nu, a,L,p)

mu = E / (2*(1+nu));
k = 3-4*nu; % plane strain
% k = (3-nu)/(1+nu); % plane stress

function [xx,yy,t1,t2] = conv_coord(x,y)
    xx = L - y;
    yy = L - x;
    t1 = atan2(yy+a, xx);
    t2 = atan2(yy-a, xx);
end

function [sxx,syy,sxy] = sol_stress(x, y)
    [x,y,t1,t2] = conv_coord(x,y);
    
    sin2t1 = sin(2*t1);
    sin2t2 = sin(2*t2);
    cos2t1 = cos(2*t1);
    cos2t2 = cos(2*t2);
    
    s_xx = p/(2*pi) .* (2*(t2-t1) + sin2t2 - sin2t1);
    s_yy = p/(2*pi) .* (2*(t2-t1) - sin2t2 + sin2t1);
    s_xy = p/(2*pi) .* (cos2t1 - cos2t2);
    
    sxx = s_yy;
    syy = s_xx;
    sxy = s_xy;
end
function [sxx] = sol_sxx(x, y)
    [sxx,~,~] = sol_stress(x,y);
end
function [syy] = sol_syy(x, y)
    [~,syy,~] = sol_stress(x,y);
end
function [sxy] = sol_sxy(x, y)
    [~,~,sxy] = sol_stress(x,y);
end

function [ux,uy] = sol_disp(x, y)
    [x,y,t1,t2] = conv_coord(x,y);
    r1 = hypot(x, y+a);
    r2 = hypot(x, y-a);
    r1(find(r1==0)) = 1.0e-8;
    r2(find(r2==0)) = 1.0e-8;
    
    u_x = (k-1).*x.*(t2-t1) + (k+1).*(y.*log(r2./r1) - a.*log(r2./(2*a)) - a.*log(r1./(2*a)));
    u_x = p/(4*pi*mu) .* u_x;
    u_y = (k-1).*((y-a).*t2 - (y+a).*t1) - (k+1).*x.*log(r2./r1);
    u_y = p/(4*pi*mu) .* u_y;
    
    ux = -u_y;
    uy = -u_x;
end
function [ux] = sol_ux(x,y)
    [ux,~] = sol_disp(x,y);
end
function [uy] = sol_uy(x,y)
    [~,uy] = sol_disp(x,y);
end



%
sol = struct();
sol.disp = @sol_disp;
sol.stress = @sol_stress;
sol.sxx = @sol_sxx;
sol.syy = @sol_syy;
sol.sxy = @sol_sxy;
sol.ux = @sol_ux;
sol.uy = @sol_uy;

return
end




