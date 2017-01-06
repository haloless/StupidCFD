function [wall] = MakeWall(xw,yw,nx,ny)

wall = struct();

wall.type = -1;

wall.xwall = [xw;yw];

nvec = [nx;ny];
wall.nwall = nvec ./ norm(nvec);

wall.dwall = dot(wall.xwall,wall.nwall);

return
end

