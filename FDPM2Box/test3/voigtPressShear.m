function [p,s] = voigtPressShear(sigma)

bm1 = [ 1; 1; 0; 1 ];

p = mean(sigma([1 2 4]));
s = sigma - p.*bm1;


return
end


