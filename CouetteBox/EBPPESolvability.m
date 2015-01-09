
function [adiag,bx,by,rhs] = EBPPESolvability(nx,ny,dx,dy, adiag,bx,by, rhs)

rhs = reshape(rhs, nx,ny);

for i = 1:nx
for j = 1:ny
    if (adiag(i,j) == 0)
        if (bx(i+1,j)==0 && bx(i,j)==0 && by(i,j+1)==0 && by(i,j)==0)
            adiag(i,j) = 1.0;
            rhs(i,j) = 0.0;
        end
    end
end
end

rhs = reshape(rhs, nx*ny,1);

return
end


