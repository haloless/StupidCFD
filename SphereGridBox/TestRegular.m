
clear all;


a = 1.0;

% np = 64
np = 160

pos = RegularPlacePoints(np);
pos = pos.';
np = size(pos,2);

dist = [];

for i = 1:np
	hi = 99999;
	for j = 1:np
		if i == j
			continue;
		end
		
		rij = pos(:,i)-pos(:,j);
		rij = sqrt(rij' * rij);
		
		hi = min(hi,rij);
	end
	dist(i,:) = [i,hi];
end

havg = mean(dist(:,2))

scale = [ 1, 1, 1 ];
% scale = [ 2, 1, 1 ];
pos(1,:) = scale(1) * pos(1,:);
pos(2,:) = scale(2) * pos(2,:);
pos(3,:) = scale(3) * pos(3,:);

figure;
PlotSphereGrid(np,pos',scale);

