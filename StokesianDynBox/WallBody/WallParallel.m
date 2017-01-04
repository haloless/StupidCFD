%% Resistance for sphere-wall in parallel direction
%%


clear all;

nmax = 1000;

a = 1;

if (0)
	% test (O'Neill 1964)
	hs = [10.0677,3.7622,1.5431,1.1276,1.0453, ...
	1.0050,1.0032,1.0018,1.0008,1.0005,1.0002];

	data = [];

	for h = hs
		[fstar,gstar] = WallParallelTrans(h/a,nmax);
		[fxt,fxr] = WallAsymptotic(h/a);
		
		data(end+1,:) = [h,fstar,gstar,fxt,fxr];
	end
end


if (0)
	% generate table
	resdat = [];
	for h = 1.001:0.001:4.000
		[fstar,gstar] = WallParallelTrans(h/a,nmax);
		
		resdat(end+1,:) = [h,fstar,gstar];
		
		if mod(h,0.1) == 0
			disp(['h=',num2str(h)]);
		end
	end
	
	dlmwrite('tmp.csv',resdat,'precision',8);
end

if (1)
	% generate table
	resdat = [];
	for h = 1.001:0.001:4.000
		[fstar,gstar] = WallParallelRot(h/a,nmax);
		
		resdat(end+1,:) = [h,fstar,gstar];
		
		if mod(h,0.1) == 0
			disp(['h=',num2str(h)]);
		end
	end
	
	dlmwrite('tmp.csv',resdat,'precision',8);
end

if (1)
	
end












