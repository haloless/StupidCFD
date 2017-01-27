function [F] = AxisymEvolverDerivForce(bridge,np,rpin,xpin)

R1 = bridge.R1;
H = bridge.H;
X2 = bridge.X2;

dh = R1*1.0e-3;

if H < dh
	% use one-sided difference
	
	[~,~,~,ene1] = AxisymEvolver(bridge, np,rpin,xpin);
	
	bridge.H = H + dh;
	bridge.X2 = X2 + dh;
	[~,~,~,ene2] = AxisymEvolver(bridge, np,rpin,xpin);
	
	bridge.H = H + dh*2;
	bridge.X2 = X2 + dh*2;
	[~,~,~,ene3] = AxisymEvolver(bridge, np,rpin,xpin);
	
	F = (-1.5*ene1 + 2.0*ene2 - 0.5*ene3) / dh;
else
	% center difference
	
	bridge.H = H + dh;
	bridge.X2 = X2 + dh;
	[~,~,~,ene2] = AxisymEvolver(bridge, np,rpin,xpin);
	
	bridge.H = H - dh;
	bridge.X2 = X2 - dh;
	[~,~,~,ene0] = AxisymEvolver(bridge, np,rpin,xpin);
	
	F = (ene2-ene0) / (dh*2);
end

bridge.H = H;
bridge.X2 = X2;

return
end

