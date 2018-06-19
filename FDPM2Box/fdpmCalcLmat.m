function [ L ] = fdpmCalcLmat(x, out_of_plane)
%fdpmCalcLmat: Calculate the tensor L = d(logb)/db
% L is reordered by [11,22,12, 33]
% if out_of_plane, L is 4x4; otherwise, L is 3x3
%
% TODO separate the tensor function part from this
%


% func = @dlgd2;
% dfunc = @ddlgd2;
func = @log;
dfunc = @recip;


Is = [ 1, 0, 0; 0, 1, 0; 0, 0, 0.5 ];

if out_of_plane
	% [xx,yy,xy,zz]
	L = zeros(4);
else
	% [xx,yy,xy]
	L = zeros(3);
end

% 
x2 = x(1:2,1:2);
[v,d] = eig(x2);
d = diag(d);

% compute f(eigenvalue)
y = func(d);
yd = dfunc(d);

% check repeated eigenvalues
repeat = 0;
diff = abs(d(1)-d(2));
maxd = max(abs(d(1)),abs(d(2)));
if maxd > 0
	diff = diff / maxd;
end
if diff < 1.0e-6
	repeat = 1;
end

if repeat
	% two duplicated eigenvalues in-plane
	L(1:3,1:3) = yd(1) * Is;
else
	% two different eigenvalues in-plane
	a = (y(1)-y(2)) / (d(1)-d(2));
	
	% two eigenprojections
	eigproj = zeros(3,2);
	% for i = 1:2
		% eigproj(1,i) = v(1,i) * v(1,i);
		% eigproj(2,i) = v(2,i) * v(2,i);
		% eigproj(3,i) = v(1,i) * v(2,i);
	% end
	eigproj(1,:) = v(1,:) .* v(1,:);
	eigproj(2,:) = v(2,:) .* v(2,:);
	eigproj(3,:) = v(1,:) .* v(2,:);
	
	e1e1 = eigproj(:,1) * eigproj(:,1)';
	e2e2 = eigproj(:,2) * eigproj(:,2)';
	% L(1:3,1:3) = a*(Is - e1e1 - e2e2) + yd(1)*e1e1 + yd(2)*e2e2;
	L(1:3,1:3) = a*Is + (yd(1)-a)*e1e1 + (yd(2)-a)*e2e2;
	
	% for i = 1:3
	% for j = 1:3
		% L(i,j) = a*(Is(i,j) - eigproj(i,1)*eigproj(j,1) - eigproj(i,2)*eigproj(j,2)) ...
		% + yd(1)*eigproj(i,1)*eigproj(j,1) + yd(2)*eigproj(i,2)*eigproj(j,2);
	% end
	% end
end

if out_of_plane
	L(4,4) = dfunc(x(3,3));
end


return
end

function [y] = recip(x)
	y = 1.0 ./ x;
end

function [y] = ddlgd2(x)
	y = 0.5 ./ x;
	return
end

function [y] = dlgd2(x)
	y = 0.5 * log(x);
	return
end

