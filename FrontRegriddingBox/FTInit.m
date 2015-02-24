
function [] = FTInit(pt0,tri0,nbr0)

FTRegridGlobals;

np = size(pt0,1);
nelem = size(tri0,1);

% copy to storage
pt(1:np,:) = pt0;
icp(1:nelem,:) = tri0;
ine(1:nelem,:) = nbr0;

% initialize link list
ffp = 1;
lfp = np;
fep = np+1;
%
ptcon(ffp:lfp-1) = ffp+1:lfp;
ptcon(lfp) = ffp;
bptcon(ffp+1:lfp) = ffp:lfp-1;
bptcon(ffp) = lfp;
%
ptcon(fep:maxpt-1) = fep+1:maxpt;
ptcon(maxpt) = fep;
bptcon(fep+1:maxpt) = fep:maxpt-1;
bptcon(fep) = maxpt;

%
ffe = 1;
lfe = nelem;
fee = nelem+1;
%
elcon(ffe:lfe-1) = ffe+1:lfe;
elcon(lfe) = ffe;
belcon(ffe+1:lfe) = ffe:lfe-1;
belcon(ffe) = lfe;
%
elcon(fee:maxel-1) = fee+1:maxel;
elcon(maxel) = fee;
belcon(fee+1:maxel) = fee:maxel-1;
belcon(fee) = maxel;

return
end

