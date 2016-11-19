
function [ tag ] = SharpIBTag(nx,ny,sdf)

tag = zeros(nx,ny);

tag_s  = -1;
tag_sf = 0;
tag_fs = 1;
tag_f  = 2;

for j = 1:ny
for i = 1:nx
    if sdf(i,j) > 0
        flag = tag_f;
        if (i>1 & sdf(i-1,j)<=0) | (i<nx & sdf(i+1,j)<=0) | (j>1 & sdf(i,j-1)<=0) | (j<ny & sdf(i,j+1)<=0)
            flag = tag_fs;
        end
    else
        flag = tag_s;
        if (i>1 & sdf(i-1,j)>0) | (i<nx & sdf(i+1,j)>0) | (j>1 & sdf(i,j-1)>0) | (j<ny & sdf(i,j+1)>0)
            flag = tag_sf;
        end
    end
    tag(i,j) = flag;
end
end


return
end


