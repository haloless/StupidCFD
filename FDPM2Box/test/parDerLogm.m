function [L] = parDerLogm(Be)
%parDerLogm: Partial Derivative of LOGM tensor function w.r.t its argument

%
[BeVec,BePr] = eig(Be);
BePr = diag(BePr);
% BePr = [BePr(1),BePr(5),BePr(9)]';

%
L = parDerGen(Be,BeVec,BePr,log(BePr),1.0 ./ BePr);

if 1
	% validate by finite difference
	Ldiff = zeros(9,9);
	
	ind = zeros(3);
	ind(1,1) = 1;
	ind(2,2) = 2;
	ind(3,3) = 3;
	ind(1,2) = 4;
	ind(2,1) = 5;
	ind(2,3) = 6;
	ind(3,2) = 7;
	ind(3,1) = 8;
	ind(1,3) = 9;
	
	db = 1.0e-6;
	
	for i = 1:3
		B = Be;
		B(i,i) = Be(i,i) + db;
		Ap = logm(B);
		
		B = Be;
		B(i,i) = Be(i,i) - db;
		Am = logm(B);
		
		dA = (Ap-Am) / (db*2);
		
		Ldiff(:,ind(i,i)) = dA([1,5,9,2,4,6,8,7,3]);
	end
	
	for i = 1:3
	for j = 1:3
	if i~=j
		% NOTE the tensor B must be symmetric for LOGM to work correctly
		
		B = Be;
		B(i,j) = Be(i,j) + db;
		B(j,i) = Be(j,i) + db;
		Ap = logm(B);
		
		B = Be;
		B(i,j) = Be(i,j) - db;
		B(j,i) = Be(j,i) - db;
		Am = logm(B);
		
		dA = (Ap-Am) / (db*2);
		dA = 0.5 * dA;
		
		Ldiff(:,ind(i,j)) = dA([1,5,9,2,4,6,8,7,3]);
	end
	end
	end
	
	disp('numer diff =')
	disp(Ldiff)
end

return
end




