function [node] = IcosahedralPoints(factor)

% basic regular icosahedron
npoint = 12;
nedge = 30;
nface = 20;
[point,edge,face] = IcosahedralRegular();

% 1. set 12 regular points
count = 12;
node = point;

% 2. refine edges
for iedge = 1:nedge
    ia = edge(1,iedge);
    ib = edge(2,iedge);
    xa = point(:,ia);
    xb = point(:,ib);
    
    for level = 1:factor-1
        count = count + 1;
        frac = level / factor;
        xfrac = (1-frac)*xa + frac*xb;
        xfrac = mapshere(xfrac);
        node(:,count) = xfrac;
    end
end

% 3. refine faces
for iface = 1:nface
    ia = face(1,iface);
    ib = face(2,iface);
    ic = face(3,iface);
    xa = point(:,ia);
    xb = point(:,ib);
    xc = point(:,ic);
    
    for l1 = 1:factor-1
    for l2 = 1:factor-1-l1
        count = count + 1;
        f1 = l1 / factor;
        f2 = l2 / factor;
        xfrac = (1-f1-f2)*xa + f1*xb + f2*xc;
        xfrac = mapshere(xfrac);
        node(:,count) = xfrac;
    end
    end
end

return
end

function [xsph] = mapshere(xin)
rin = norm(xin);
xsph = xin ./ rin;
return
end


