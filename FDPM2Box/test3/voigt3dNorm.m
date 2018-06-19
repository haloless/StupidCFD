function [n] = voigt3dNorm(v)
%voigt3dNorm

% assert(numel(v) == 6);

n = sum(v(1:3).^2) + 2*sum(v(4:6).^2);
n = sqrt(n);


return
end


