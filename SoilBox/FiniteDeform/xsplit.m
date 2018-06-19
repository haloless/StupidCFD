function [ex,ey,ez] = xsplit(edof,coord,dof,nen,nel)
% 


nn = zeros(1,nen);

ex = zeros(nel,nen);
ey = ex;
ez = ex;

for i = 1:nel
	for j = 1:nen
		nn(j) = find(dof(:,1)-edof(i,j*3-1)==0);
	end
	ex(i,:) = coord(nn,1)';
	ey(i,:) = coord(nn,2)';
	ez(i,:) = coord(nn,3)';
end


return
end