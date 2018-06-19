function [sol] = refsol_TimoshenkoCantilever(E0,nu0, L0,c0,P0)
%refsol_TimoshenkoCantilever
% Timoshenko's cantilever
% we assume plain-strain state (which is different from plain-stress commonly referred)
%
% Inputs:
% L0: length
% c0: half width 
% P0: total load at free end
% 


% modified for plain-strain
nubar = nu0 / (1-nu0);
Ebar = E0 / (1-nu0^2);
% Ebar = E0 / (1-nubar^2);

% for plain-stress
% nubar = nu0;
% Ebar = E0;

% cantilever's sectional moment
d0 = c0 * 2;
I0 = d0^3 / 12;

%
timo_ux = @(x,y) -P0/(6*Ebar*I0).* y .*((6*L0-3*x).*x + (2+nubar).*(y.^2-c0^2));
timo_uy = @(x,y) P0/(6*Ebar*I0).*(3*nubar*y.^2 .*(L0-x) + (4+5*nubar)*c0^2*x + (3*L0-x).*x.^2);
timo_sxx = @(x,y) -P0/I0 * (L0-x).*y;
timo_syy = @(x,y) zeros(size(x));
timo_sxy = @(x,y) P0/(2*I0) * (c0^2 - y.^2);

% return an enclosure
sol = struct();
sol.ux = timo_ux;
sol.uy = timo_uy;
sol.sxx = timo_sxx;
sol.syy = timo_syy;
sol.sxy = timo_sxy;




return
end


