%% Resistance for sphere-wall in perpendicular direction
%%
%%

clear all;

nmax = 100;

a = 1;

if (0)
	% test [Brenner,1961]
	
	hs = [1.1276260,1.5430806,2.3524096,3.7621957,6.1822895,10.067662];
	data = [];
	for h = hs
		f = WallPerpendicularTrans(h/a,nmax);
		data(end+1,:) = [h,f];
	end
end

if (1)
	% generate table
	
	resdat = [];
	for h = 1.001:0.001:4.000
		f = WallPerpendicularTrans(h/a,nmax);
		g = WallPerpendicularRot(h/a,nmax);
		
		resdat(end+1,:) = [h,f,g];
		
		if mod(h,0.1) == 0
			disp(['h=',num2str(h)]);
		end
	end
	
	dlmwrite('tmp.csv',resdat,'precision',8);
end


