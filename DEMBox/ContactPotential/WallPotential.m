function [phi] = WallPotential(wall,x,y)
pos = [x';y'];
phi = wall.nwall' * pos;
phi = phi' - wall.dwall;
return
end


