
function [ kenew ] = fdelel(ke,n1)
% This function will delete KE and KEC.
% Then the input KE index will be no longer valid after that.
% So a new index KE must be returned instead. 

FTRegridGlobals;


kenew = ke;

[n1,n2,n3,nc1,nc2,nc3,kp] = flocal(ke,n1,icp,ine);
kec = ine(ke,n1);

% check if the delete operation will result in
% bad mesh configuration (@see fdelet.F)
bad_delete = 0;
ke2 = ine(ke,n2);
ke3 = ine(ke,n3);
kec1 = ine(kec,nc1);
kec2 = ine(kec,nc2);
for i = 1:3
    el1 = ine(kec1,i);
    el2 = ine(ke2,i);
    for j = 1:3
        if el1==ine(kec2,j)
            if (el1~=kec); bad_delete=1; break; end
        end
        if el2==ine(ke3,j)
            if (el2~=ke); bad_delete=1; break; end
        end
    end
    if bad_delete; break; end
end
% TODO check kp(5:8)

if bad_delete; return; end

% set p1 position
pt(kp(1),:) = 0.5*(pt(kp(1),:)+pt(kp(2),:));
% eliminate p2
if kp(2)==ffp
    ffp = ptcon(kp(2));
elseif kp(2)==lfp
    lfp = bptcon(kp(2));
end
np = np - 1;
ptcon(bptcon(kp(2))) = ptcon(kp(2));
bptcon(ptcon(kp(2))) = bptcon(kp(2));
bptcon(fep) = kp(2);
ptcon(kp(2)) = fep;
fep = kp(2);

%
kenew = elcon(ke);
if kenew == kec
    kenew = elcon(kec);
end

% eliminate KE and KEC
[nelem,ffe,lfe,fee,elcon,belcon] = LinkListDel(ke, ...
nelem,ffe,lfe,fee,elcon,belcon, maxel);
[nelem,ffe,lfe,fee,elcon,belcon] = LinkListDel(kec, ...
nelem,ffe,lfe,fee,elcon,belcon, maxel);

% set element pointers
for i = 1:3
    if ine(ine(ke,n2),i) == ke
        ine(ine(ke,n2),i) = ine(ke,n3);
    end
    if ine(ine(ke,n3),i) == ke
        ine(ine(ke,n3),i) = ine(ke,n2);
    end
    if ine(ine(kec,nc2),i) == kec
        ine(ine(kec,nc2),i) = ine(kec,nc1);
    end
    if ine(ine(kec,nc1),i) == kec
        ine(ine(kec,nc1),i) = ine(kec,nc2);
    end
end

% redirect pt(2) to pt(1)
% this requires a loop through all elements holding pt(2) previously
ke1 = ine(ke,n2);
while (1)
    mo = 0;
    for i = 1:3
        if icp(ke1,i) == kp(2)
            icp(ke1,i) = kp(1);
            mo = i;
        end
    end
    
    if ke1 == ine(kec,nc2); break; end
    ke1 = ine(ke1,mo);
end

return
end


