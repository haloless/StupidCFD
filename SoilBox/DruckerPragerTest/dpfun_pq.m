function [p,q] = dpfun_pq(sigma)
% SIGMA must be principle stress
% p is hydrostatic mean stress
% q is sqrt(3/2) * |s|

s1 = sigma(1);
s2 = sigma(2);
s3 = sigma(3);

p = (s1+s2+s3) / 3;
q = sqrt(s1^2 + s2^2 + s3^2 - s1*s2 - s2*s3 - s3*s1);


return
end

