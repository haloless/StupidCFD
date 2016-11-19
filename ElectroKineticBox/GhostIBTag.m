
function [ tag ] = GhostIBTag(nx,ny,sdf)

tag = zeros(nx,ny);
tag(sdf>0) = 1;
tag(sdf<=0) = -1;
for j = 2:ny-1
for i = 2:nx-1
    if sdf(i,j) <= 0
        if sdf(i-1,j)>0 | sdf(i+1,j)>0 | sdf(i,j-1)>0 | sdf(i,j+1)>0
            tag(i,j) = 0;
        end
    end
end
end

return
end


