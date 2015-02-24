
function [] = faddel(ke,n1)

FTRegridGlobals;



[n1,n2,n3,nc1,nc2,nc3,kp] = flocal(ke,n1,icp,ine);

kec = ine(ke,n1);

% add new point and put in list
if (0) 
    newpt = fep;
    fep = ptcon(fep);
    ptcon(newpt) = ptcon(lfp);
    ptcon(lfp) = newpt;
    lfp = newpt;
    np = np + 1;
else
    [newpt,np,ffp,lfp,fep,ptcon,bptcon] = LinkListAdd( ...
    np,ffp,lfp,fep,ptcon,bptcon,maxpt);
end
% set point position
% TODO use a high-order interpolation
pt(newpt,:) = 0.5 * (pt(kp(1),:)+pt(kp(2),:));

% add one element and put in list
if (0)
    newel = fee;
    fee = elcon(fee);
    elcon(newel) = elcon(lfe);
    elcon(lfe) = newel;
    lfe = newel;
    nelem = nelem + 1;
else
    [newel,nelem,ffe,lfe,fee,elcon,belcon] = LinkListAdd( ...
    nelem,ffe,lfe,fee,elcon,belcon,maxel);
end
% add another
if (0)
    newelc = fee;
    fee = elcon(fee);
    elcon(newelc) = elcon(lfe);
    elcon(lfe) = newelc;
    lfe = newelc;
    nelem = nelem + 1;
else
    [newelc,nelem,ffe,lfe,fee,elcon,belcon] = LinkListAdd( ...
    nelem,ffe,lfe,fee,elcon,belcon,maxel);
end

% modify current element
icp(ke,n1) = newpt;
% modify current neighbor
icp(kec,nc1) = newpt;
% set new element
icp(newel,1) = newpt;
icp(newel,2) = kp(3);
icp(newel,3) = kp(1);
% set new neighbor
icp(newelc,1) = newpt;
icp(newelc,2) = kp(1);
icp(newelc,3) = kp(4);

% adjust neighborhood
ke3 = ine(ke,n3);
kec1 = ine(kec,nc1);
ine(ke,n3) = newel;
ine(kec,nc1) = newelc;
% set neighborhood
ine(newel,1) = ke;
ine(newel,2) = ke3;
ine(newel,3) = newelc;
ine(newelc,1) = newel;
ine(newelc,2) = kec1;
ine(newelc,3) = kec;

% don't forget these two neighbor...
for i = 1:3
    if ine(ke3,i)==ke
        ine(ke3,i) = newel;
    end
    if ine(kec1,i)==kec
        ine(kec1,i) = newelc;
    end
end


return
end

