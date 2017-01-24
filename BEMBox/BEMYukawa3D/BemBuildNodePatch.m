% build by Node Patch

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
		
		for ind = 1:3
			j = face(jelem,ind);
			if ind==1
				j1 = face(jelem,2);
				j2 = face(jelem,3);
			elseif ind==2
				j1 = face(jelem,3);
				j2 = face(jelem,1);
			else
				j1 = face(jelem,1);
				j2 = face(jelem,2);
			end
			
			level = 1;
			if i == j
				level = 3;
			end
			
			vj = vert(:,j);
			va = (vert(:,j)+vert(:,j1))/2;
			vb = (vert(:,j)+vert(:,j2))/2;
			
			vp = vert(:,i);
			% vp = facecent(:,i);
			
			[aa,bb] = BemCoefNodePatch(vp, ...
				jelem, j, vj,va,vb, level);
			if i == j
				bb = bb + 0.5/patchcnt(i);
			end
			if bc(1,j) == 1
				amat(i,j) = amat(i,j) - aa;
				bvec(i) = bvec(i) - bb*bc(2,j);
			elseif bc(1,j) == 2
				amat(i,j) = amat(i,j) + bb;
				bvec(i) = bvec(i) + aa*bc(2,j);
			end
		end
		
		% level = 1;
		
		% % node 1
		% [aa,bb] = BemCoefNodePatch(facecent(:,i), ...
			% jelem, j1, v1, v4, v6, level);
		% if i == j1
			% bb = bb + 0.5/patchcnt(i);
		% end
		% if bc(1,j1) == 1
			% amat(i,j1) = amat(i,j1) - aa;
			% bvec(i) = bvec(i) - bb*bc(2,j1);
		% elseif bc(1,j1) == 2
			% amat(i,j1) = amat(i,j1) + bb;
			% bvec(i) = bvec(i) + aa*bc(2,j1);
		% end
		
		% % node 2
		% [aa,bb] = BemCoefNodePatch(facecent(:,i), ...
			% jelem, j2, v2, v5, v4, level);
		% if i == j2
			% bb = bb + 0.5/patchcnt(i);
		% end
		% if bc(1,j2) == 1
			% amat(i,j2) = amat(i,j2) - aa;
			% bvec(i) = bvec(i) - bb*bc(2,j2);
		% elseif bc(1,j2) == 2
			% amat(i,j2) = amat(i,j2) + bb;
			% bvec(i) = bvec(i) + aa*bc(2,j2);
		% end
		
		% % node 3
		% [aa,bb] = BemCoefNodePatch(facecent(:,i), ...
			% jelem, j3, v3, v6, v5, level);
		% if i == j3
			% bb = bb + 0.5/patchcnt(i);
		% end
		% if bc(1,j3) == 1
			% amat(i,j3) = amat(i,j3) - aa;
			% bvec(i) = bvec(i) - bb*bc(2,j3);
		% elseif bc(1,j3) == 2
			% amat(i,j3) = amat(i,j3) + bb;
			% bvec(i) = bvec(i) + aa*bc(2,j3);
		% end
		
		
		% % check distance for quadrature rule
		% rij = norm(facecent(:,i)-facecent(:,j));
		% dlen = facelmax(i) + facelmax(j);
		% level = 0;
		% if i==j || rij<dlen*1.0e-3
			% % singular integral
			% level = 3;
		% elseif rij < dlen*2
			% % close integral
			% level = 2;
		% else
			% level = 1;
		% end
		
		% [aa,bb] = BemCoef2(facecent(:,i), j, level);
		% if i == j
			% bb = bb + 0.5;
		% end
		
		% if bc(1,j) == 1
			% amat(i,j) = amat(i,j) - aa;
			% bvec(i) = bvec(i) - bb*bc(2,j);
		% elseif bc(1,j) == 2
			% amat(i,j) = amat(i,j) + bb;
			% bvec(i) = bvec(i) + aa*bc(2,j);
		% end
	end
	if mod(jelem,20) == 0
		disp(['elem=',int2str(jelem)]);
	end
end


