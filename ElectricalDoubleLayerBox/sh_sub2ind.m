function [ind] = sh_sub2ind(n,m, islocal)
%

ind = m + n + 1;

if nargin == 2
	ind = ind + n^2;
end

return
end
