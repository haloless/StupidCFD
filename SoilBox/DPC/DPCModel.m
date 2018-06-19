function [p,s] = DPCModel(p,s)

DPCGlobals;

ltype = 1;
if layer >= istat
    ltype = 2;
end

% intersection 
fcut = log(aa/ac) / ab;
% tension cutoff
tcut = 0.0;

sigx = s(1) + p;
sigy = s(2) + p;
sigz = 3*p - sigx - sigy;
sigxz = s(3);
sigyz = s(4);
sigxy = s(5);
% 
ep = el;

DPCMatLaw;

%
p = (sigx+sigy+sigz) / 3;
s(1) = sigx - p;
s(2) = sigy - p;
s(3) = sigxz;
s(4) = sigyz;
s(5) = sigxy;
el = ep;



return
end


