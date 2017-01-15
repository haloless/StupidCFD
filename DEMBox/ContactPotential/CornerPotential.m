function [phi] = CornerPotential(corner,x,y)

pos = [x;y];

vw = pos - corner.xw;

dot1 = dot(vw, corner.n1);
dot2 = dot(vw, corner.n2);

proj1 = pos - dot1*corner.n1;
proj2 = pos - dot2*corner.n2;

para1 = proj1 - corner.xw;
para2 = proj2 - corner.xw;

dist = 9999999;
if dot(para1,corner.t1) < 0
    if abs(dist) > abs(dot1)
        dist = dot1;
    end
end
if dot(para2,corner.t2) > 0
    if abs(dist) > abs(dot2)
        dist = dot2;
    end
end
if 1
    distw = norm(vw);
    if abs(dist) > distw
        signw = sign(dot(vw,corner.nw));
        dist = distw * signw;
    end
end

phi = dist;

return
end



