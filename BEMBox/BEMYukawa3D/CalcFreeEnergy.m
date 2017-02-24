
%
% Free energy for constant potential case
% TODO constant charge case
%

ene = 0;
for i = 1:nface
	i1 = face(i,1);
	i2 = face(i,2);
	i3 = face(i,3);
	
	un = (sol(i1)+sol(i2)+sol(i3)) / 3.0;
	u = bc(2,i);
	
	ds = facearea(i);
	
	ene = ene + un*u*ds;
end

ene = 0.5 * ene;

