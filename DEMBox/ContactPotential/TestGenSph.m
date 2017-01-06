
clear all;

% mdiv = 8;
% ndiv = 6;
mdiv = 12;
ndiv = 8;

dphi = 2*pi / mdiv;
dtheta = pi / ndiv;

if 0
nf = mdiv*ndiv;
np = mdiv*(ndiv-1) + 2;

xp = zeros(np,3);
xp(1,:) = [0,0,1];
xp(np,:) = [0,0,-1];
cnt = 2;
for m = 1:mdiv
	phi = (m-1)*dphi;
	sphi = sin(phi);
	cphi = cos(phi);
	
	for n = 1:ndiv-1
		theta = n*dtheta;
		st = sin(theta);
		ct = cos(theta);
		
		xp(cnt,:) = [cphi*st,sphi*st,ct];
		cnt = cnt + 1;
	end
end

el = zeros(nf,3);
for m = 1:mdiv
for n = 1:ndiv
	
end
end

end

nf = mdiv * ndiv;
np = mdiv * (ndiv+1);

xp = zeros(np,3);
for m = 1:mdiv
	phi = (m-1)*dphi;
	sphi = sin(phi);
	cphi = cos(phi);
	for n = 1:ndiv+1
		theta = (n-1)*dtheta;
		st = sin(theta);
		ct = cos(theta);
		
		idx = (m-1)*(ndiv+1) + n;
		xp(idx,:) = [cphi*st,sphi*st,ct];
	end
end

el = zeros(nf,4);
for m = 1:mdiv
	for n = 1:ndiv
		i1 = (m-1)*(ndiv+1) + n;
		i2 = i1 + 1;
		i3 = (mod(m,mdiv))*(ndiv+1) + n;
		i4 = i3 + 1;
		
		idx = (m-1)*ndiv + n;
		el(idx,:) = [i1,i2,i4,i3];
	end
end

fileid = fopen('tmp.obj','w');
fprintf(fileid, '#test\n');

for i = 1:np
	fprintf(fileid, 'v %f %f %f\n', xp(i,1),xp(i,2),xp(i,3));
end

for i = 1:nf
	fprintf(fileid, 'f %d %d %d %d\n', el(i,1),el(i,2),el(i,3),el(i,4));
end


fclose(fileid);

% figure;
% plot3(xp(:,1),xp(:,2),xp(:,3),'x');
% axis equal;

