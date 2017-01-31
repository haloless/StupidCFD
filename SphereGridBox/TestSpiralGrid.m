
clear all;


a = 1.0;

% h = a / 4.0
h = a * 0.32

np = SpiralPointsEstimNum(h/a)

np = 64

pos = SpiralPoints(a,[0,0,0]',np);

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

scale = [ 2, 1, 1 ];
pos(1,:) = scale(1) * pos(1,:);
pos(2,:) = scale(2) * pos(2,:);
pos(3,:) = scale(3) * pos(3,:);

figure;
PlotSphereGrid(np,pos',scale);

