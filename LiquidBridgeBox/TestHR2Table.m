%%
%% 
%%

clear all;

sigma = 1.0;
R = 1.0;
theta = deg2rad(0);
R1 = R;
% R2 = R;
R2 = -1;
theta1 = theta;
theta2 = theta;

vol = [ ...
1.0e-7, 2.0e-7, 3.0e-7, 4.0e-7, 5.0e-7, 6.0e-7, 7.0e-7, 8.0e-7, 9.0e-7, ...
1.0e-6, 2.0e-6, 3.0e-6, 4.0e-6, 5.0e-6, 6.0e-6, 7.0e-6, 8.0e-6, 9.0e-6, ...
1.0e-5, 2.0e-5, 3.0e-5, 4.0e-5, 5.0e-5, 6.0e-5, 7.0e-5, 8.0e-5, 9.0e-5, ...
1.0e-4, 2.0e-4, 3.0e-4, 4.0e-4, 5.0e-4, 6.0e-4, 7.0e-4, 8.0e-4, 9.0e-4, ...
1.0e-3, 2.0e-3, 3.0e-3, 4.0e-3, 5.0e-3, 6.0e-3, 7.0e-3, 8.0e-3, 9.0e-3, ...
1.0e-2, 2.0e-2, 3.0e-2, 4.0e-2, 5.0e-2, 6.0e-2, 7.0e-2, 8.0e-2, 9.0e-2, ...
1.0e-1, 1.5e-1, 2.0e-1, 2.5e-1, 3.0e-1, ...
];

% regularized by rupture distance
ndist = 101;
dist = linspace(0,1,ndist);

data = [];
for V = vol
	if R2 > 0
		hrup = BridgeRuptureLian(theta,V);
	else
		hrup = BridgeRuptureMikami(R1,R2,theta,V);
	end
	
	data1 = [];
	for hstar = dist
		% recover true distance
		H = hrup * hstar;
		
		data1(end+1) = BridgeForceHR2(R1,R2,H,theta1,theta2,V,sigma);
	end
	
	data0 = [];
	if 0
		data0 = zeros(size(data1));
		for idist = 1:ndist
			hstar = dist(idist);
			H = hrup * hstar;
			
			if mod(hstar,0.2) ~= 0
				continue;
			end
			
			np = 101;
			bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
			data0(idist) = AxisymEvolverDirectForce(bridge,np,[],[]);
		end
	end
	
	data = [data, data0',data1'];
end

figure;
plot(data,'x-');
axis([0 60 0 7]);






