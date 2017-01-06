
function [mat] = WallResistanceMatrix(z)

global wallresdat;

if isempty(wallresdat)
	% h/a,Fxt,Fxr,Tyt,Tyr,Fzt,Tzr
	wallresdat = dlmread('WallSphereResTab.csv');
	disp(['Load tabular data: ', int2str(size(wallresdat))]);
end

zvalid = 4.0;
zsmall = 1.001;

if z <= zsmall
	% asymptotic lubrication theory
	[fxt,fxr,tyr,fzt,tzr] = WallAsymptotic(z);
	tyt = fxr;
elseif z <= zvalid
	resint = interp1(wallresdat(:,1),wallresdat(:,2:end),z);
	fxt = resint(1);
	fxr = resint(2);
	tyt = resint(3);
	tyr = resint(4);
	fzt = resint(5);
	tzr = resint(6);
else
	error(['Invalid sphere-wall distance z=',num2str(z)]);
end

mat = [...
fxt, 0, 0, 0, -fxr, 0;
0, fxt, 0, fxr, 0, 0;
0, 0, fzt, 0, 0, 0;
0, tyt, 0, tyr, 0, 0;
-tyt, 0, 0, 0, tyr, 0;
0, 0, 0, 0, 0, tzr];


return
end
