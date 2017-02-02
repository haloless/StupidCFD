% build by linear element

patchcnt = zeros(nvert,1);
for jelem = 1:nface
	% face nodes
	j1 = face(jelem,1);
	j2 = face(jelem,2);
	j3 = face(jelem,3);
	
	patchcnt(j1) = patchcnt(j1) + 1;
	patchcnt(j2) = patchcnt(j2) + 1;
	patchcnt(j3) = patchcnt(j3) + 1;
end

amat = zeros(nvert,nvert);
bvec = zeros(nvert,1);
for jelem = 1:nface
	
	% face nodes
	j1 = face(jelem,1);
	j2 = face(jelem,2);
	j3 = face(jelem,3);
	
	% node position
	v1 = vert(:,j1);
	v2 = vert(:,j2);
	v3 = vert(:,j3);
	v4 = (v1+v2)/2;
	v5 = (v2+v3)/2;
	v6 = (v3+v1)/2;
	
	
	% face center
	vcentj = facecent(:,jelem);
	% face normal
	vnormj = facenvec(:,jelem);
	
	for i = 1:nvert
		vp = vert(:,i);
        
        level = 1;
        j1 = face(jelem,1);
        j2 = face(jelem,2);
        j3 = face(jelem,3);
        if face(jelem,1) == i
            level = 3;
            j1 = face(jelem,1);
            j2 = face(jelem,2);
            j3 = face(jelem,3);
        elseif face(jelem,2) == i
            level = 3;
            j1 = face(jelem,2);
            j2 = face(jelem,3);
            j3 = face(jelem,1);
        elseif face(jelem,3) == i
            level = 3;
            j1 = face(jelem,3);
            j2 = face(jelem,1);
            j3 = face(jelem,2);
        end
        
        v1 = vert(:,j1);
        v2 = vert(:,j2);
        v3 = vert(:,j3);
        
        [aa,bb] = BemCoefLinear(vp,jelem,v1,v2,v3, level);
        
        if i == j1
            bb(1) = bb(1) + 0.5/patchcnt(i);
        end
        if bc(1,jelem) == 1
            amat(i,j1) = amat(i,j1) - aa(1);
            bvec(i) = bvec(i) - bb(1)*bc(2,j1);
            amat(i,j2) = amat(i,j2) - aa(2);
            bvec(i) = bvec(i) - bb(2)*bc(2,j2);
            amat(i,j3) = amat(i,j3) - aa(3);
            bvec(i) = bvec(i) - bb(3)*bc(2,j3);
            
        elseif bc(1,jelem) == 2
            % amat(i,j) = amat(i,j) + bb;
            % bvec(i) = bvec(i) + aa*bc(2,j);
            error('not supported');
        end
        
	end
	if mod(jelem,20) == 0
		disp(['elem=',int2str(jelem)]);
	end
end


