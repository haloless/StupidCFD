
fbnd = zeros(numDofs,1);
kbnd = zeros(numDofs,1);

% nodeCurr = nodeCoord + reshape(uv, 2,numNodes);

for ibndry = 1:nbndry
    p1 = bndry(ibndry).p1;
    p2 = bndry(ibndry).p2;
    tlen = norm(p2-p1);
    tvec = bndry(ibndry).tvec;
    nvec = bndry(ibndry).nvec;
    aa = bndry(ibndry).a;
    bb = bndry(ibndry).b;
    cc = bndry(ibndry).c;
    
    for i = 1:numNodes
        ipos = nodeCurr(:,i);
        tproj = dot(ipos-p1,tvec) / tlen;
        ndist = dot(ipos-p1,nvec);
        
        rb = h0/2;
        if 0<=tproj && tproj<=1 && abs(ndist)<rb*2
            % i
            ovlp = rb - ndist; % assert(ovlp >= 0);
            
            idof = fdpmNodeDof(i, 'xy');
            
            fbnd(idof) = fbnd(idof) + kbndry*ovlp*nvec;
            
            kbnd(idof) = kbnd(idof) + kbndry*nvec.^2;
        end
    end
end


