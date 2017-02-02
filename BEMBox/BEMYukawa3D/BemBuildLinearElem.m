% build by linear element


% count elements for each node
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

disp('dihedral angle');
%
dihedral_angle = zeros(nedge,1);
for iedge = 1:nedge
	% two nodes
	v1 = vert(:,edge(iedge,1));
	v2 = vert(:,edge(iedge,2));
	
	% two elements
	ja = edge2face(iedge,1);
	jb = edge2face(iedge,2);
	
	% normal vector
	na = facenvec(:,ja);
	nb = facenvec(:,jb);
	
	beta = pi - acos(dot(na,nb));
	
	% edge normal
	% ne = na + nb; ne = ne ./ norm(ne);
	% ve = (v1+v2) / 2;
	
	% correct to be internal dihedral angle
	te = v2 - v1;
	ta = cross(na,te);
	tb = cross(te,nb);
	tt = cross(ta,tb);
	check = dot(tt,te);
	if check > 0
		beta = pi*2 - beta;
	end
	
	dihedral_angle(iedge) = beta;
end

disp('solid angle');
% calculate solid angle for node
solid_angle = zeros(nvert,1);
solid_angle(:) = pi*2;
for iedge = 1:nedge
	% two nodes
	i1 = edge(iedge,1);
	i2 = edge(iedge,2);
	
	beta = dihedral_angle(iedge);
	solid_angle(i1) = solid_angle(i1) + beta - pi;
	solid_angle(i2) = solid_angle(i2) + beta - pi;
end
% convert to external solid angle and regularize by 4pi
for ivert = 1:nvert
	ang = solid_angle(ivert);
	solid_angle(ivert) = 1 - ang/(pi*4);
end


amat = zeros(nvert,nvert);
bvec = zeros(nvert,1);
for jelem = 1:nface
	
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
        
		if level == 99
			% fprintf('vert=%d,elem=%d,j1=%d,j2=%d,j3=%d\n',i,jelem,j1,j2,j3);
		end
		
        [aa,bb] = BemCoefLinear(vp, jelem,v1,v2,v3, level);
        
        if i == j1
			% sangle = 0.5;
			% sangle = 0.55;
			sangle = solid_angle(i);
            bb(1) = bb(1) + sangle/patchcnt(i);
        end
        if bc(1,jelem) == 1
            amat(i,j1) = amat(i,j1) - aa(1);
            bvec(i) = bvec(i) - bb(1)*bc(2,j1);
			
            amat(i,j2) = amat(i,j2) - aa(2);
            bvec(i) = bvec(i) - bb(2)*bc(2,j2);
			
            amat(i,j3) = amat(i,j3) - aa(3);
            bvec(i) = bvec(i) - bb(3)*bc(2,j3);
            
        elseif bc(1,jelem) == 2
			amat(i,j1) = amat(i,j1) + bb(1);
			bvec(i) = bvec(i) + aa(1)*bc(2,j1);
			amat(i,j2) = amat(i,j2) + bb(2);
			bvec(i) = bvec(i) + aa(2)*bc(2,j2);
			amat(i,j3) = amat(i,j3) + bb(3);
			bvec(i) = bvec(i) + aa(3)*bc(2,j3);
			
            % amat(i,j) = amat(i,j) + bb;
            % bvec(i) = bvec(i) + aa*bc(2,j);
            % error('not supported');
        end
        
	end
	if mod(jelem,20) == 0
		disp(['elem=',int2str(jelem)]);
	end
end


