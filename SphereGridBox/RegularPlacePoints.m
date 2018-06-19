
function [xg] = RegularPlacePoints(ng)


r = 1.0;
a = 4*pi*r^2 / ng;
d = sqrt(a);

mt = round(pi/d);
dt = pi / mt;
dp = a / dt;

% xg = zeros(ng, 3);
xg = [];
cnt = 0;
for m = 0:(mt-1)
    t = pi * (m+0.5) / mt;
    mp = round(2*pi*sin(t)/dp);
    
    for n = 0:(mp-1)
        p = 2*pi*n/mp;
        x = r * sin(t) * cos(p);
        y = r * sin(t) * sin(p);
        z = r * cos(t);
        xg(end+1,:) = [x,y,z];
        cnt = cnt + 1;
    end
end



return
end


