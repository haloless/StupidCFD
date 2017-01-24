% build by Constant Element 

amat = zeros(nface,nface);
bvec = zeros(nface,1);
for j = 1:nface
	for i = 1:nface
		% check distance for quadrature rule
		rij = norm(facecent(:,i)-facecent(:,j));
		dlen = facelmax(i) + facelmax(j);
		level = 0;
		if i==j || rij<dlen*1.0e-3
			% singular integral
			level = 3;
		elseif rij < dlen*2
			% close integral
			level = 2;
		else
			level = 1;
		end
		
		[aa,bb] = BemCoef2(facecent(:,i), j, level);
		if i == j
			bb = bb + 0.5;
		end
		
		if bc(1,j) == 1
			amat(i,j) = amat(i,j) - aa;
			bvec(i) = bvec(i) - bb*bc(2,j);
		elseif bc(1,j) == 2
			amat(i,j) = amat(i,j) + bb;
			bvec(i) = bvec(i) + aa*bc(2,j);
		end
	end
	if mod(j,20) == 0
		disp(['elem=',int2str(j)]);
	end
end


