
function [ A, rb ] = LaplacianOp2D(nx,ny,dx,dy, bctype,bcvalue)

bc_neu = 0;
bc_dir = 1;
bc_per = 2;

bc_xlo = bctype(1,1);
bc_xhi = bctype(1,2);
bc_ylo = bctype(2,1);
bc_yhi = bctype(2,2);
uxlo = bcvalue(1,1);
uxhi = bcvalue(1,2);
uylo = bcvalue(2,1);
uyhi = bcvalue(2,2);


dx2 = dx^2;
dy2 = dy^2;
rdx2 = 1.0 / dx2;
rdy2 = 1.0 / dy2;

nall = nx*ny;
A = sparse(nall,nall);
rb = zeros(nall,1);

ind = zeros(nx,ny);
ind(:) = 1:nall;

for j = 1:ny
for i = 1:nx
    idx = ind(i,j);
    
    aa = 0;
    
    if i > 1
        cc = rdx2;
        aa = aa - cc;
        A(idx,ind(i-1,j)) = cc;
    else
        if bc_xlo == bc_neu
            rb(idx) = rb(idx) - uxlo/dx;
        elseif bc_xlo == bc_dir
            cc = 2.0*rdx2;
            aa = aa - cc;
            rb(idx) = rb(idx) + uxlo*cc;
        elseif bc_xlo == bc_per
            cc = rdx2;
            aa = aa - cc;
            A(idx,ind(nx,j)) = cc;
        end
    end
    
    if i < nx
        cc = rdx2;
        aa = aa - cc;
        A(idx,ind(i+1,j)) = cc;
    else
        if bc_xhi == bc_neu
            rb(idx) = rb(idx) + uxhi/dx;
        elseif bc_xhi == bc_dir
            cc = 2.0*rdx2;
            aa = aa - cc;
            rb(idx) = rb(idx) + uxhi*cc;
        elseif bc_xhi == bc_per
            cc = rdx2;
            aa = aa - cc;
            A(idx,ind(1,j)) = cc;
        end
    end
    
    if j > 1
        cc = rdy2;
        aa = aa - cc;
        A(idx,ind(i,j-1)) = cc;
    else
        if bc_ylo == bc_neu
            rb(idx) = rb(idx) - uylo/dy;
        elseif bc_ylo == bc_dir
            cc = 2.0*rdy2;
            aa = aa - cc;
            rb(idx) = rb(idx) + uylo*cc;
        elseif bc_ylo == bc_per
            cc = rdx2;
            aa = aa - cc;
            A(idx,ind(i,ny)) = cc;
        end
    end
    
    if j < ny
        cc = rdy2;
        aa = aa - cc;
        A(idx,ind(i,j+1)) = cc;
    else
        if bc_yhi == bc_neu
            rb(idx) = rb(idx) + uyhi/dy;
        elseif bc_yhi == bc_dir
            cc = 2.0*rdy2;
            aa = aa - cc;
            rb(idx) = rb(idx) + uyhi*cc;
        elseif bc_yhi == bc_per
            cc = rdy2;
            aa = aa - cc;
            A(idx,ind(i,1)) = cc;
        end
    end
    
    A(idx,idx) = aa;
end
end



return
end


