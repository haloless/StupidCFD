function [K] = assem(edof,K,Ke)
% assemble element stiffness to global stiffness

[nie,n] = size(edof);

t = edof(:,2:n);

for i = 1:nie
	ind = t(i,:);
	K(ind,ind) = K(ind,ind) + Ke;
end



return
end


