function [ia,ib,ic, la,lb,lc, ta,tb,tc] = triRotateMaxAngle(tri, pos, i)
%triRotateMaxAngle: Rotate index

i1 = tri(i,1); 
i2 = tri(i,2); 
i3 = tri(i,3);

l1 = norm(pos(i2,:)-pos(i3,:));
l2 = norm(pos(i3,:)-pos(i1,:));
l3 = norm(pos(i1,:)-pos(i2,:));

ll = [l1,l2,l3];

ii = 1;
if l2 > ll(ii)
    ii = 2;
end
if l3 > ll(ii)
    ii = 3;
end

if ii == 1
    ia = i1; ib = i2; ic = i3;
    la = l1; lb = l2; lc = l3;
elseif ii == 2
    ia = i2; ib = i3; ic = i1;
    la = l2; lb = l3; lc = l1;
else
    ia = i3; ib = i1; ic = i2;
    la = l3; lb = l1; lc = l2;
end
assert(la>=lb && la>=lc);

if nargout > 6
    ta = acos((lb^2+lc^2-la^2)/(2*lb*lc));
    tb = acos((lc^2+la^2-lb^2)/(2*lc*la));
    tc = pi - ta - tb;
end




return
end

