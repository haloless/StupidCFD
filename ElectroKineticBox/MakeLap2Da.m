
function [ Lap, rb ] = MakeLap2Da(nx,ny,dx,dy, bctype,bcval)
% Return matrix and BC part, so that L(x) = A.x + r
%

bc_neu = 0;
bc_dir = 1;
bc_per = 2;

% check periodic BC match
if bctype(1,1)==bc_per | bctype(1,2)==bc_per
    if bctype(1,1) ~= bctype(1,2)
        error('PBC in x not match');
    end
end
if bctype(2,1)==bc_per | bctype(2,2)==bc_per
    if bctype(2,1) ~= bctype(2,2)
        error('PBC in y not match');
    end
end


rdx = 1/dx;
rdy = 1/dy;
rdx2 = 1/dx^2;
rdy2 = 1/dy^2;

np = nx * ny;
nbuf = np * 5;

ind = reshape(1:np, nx,ny);
tuples = [];
rb = zeros(np,1);

%
% the internal bulk matrix
%
I = 2:nx-1;
J = 2:ny-1;
val = zeros(nx-2,ny-2);

idp = ind(I,J);
idn = ind(I-1,J);
val(:) = rdx2;
tuples = [tuples; idp(:), idn(:), val(:) ];

idp = ind(I,J);
idn = ind(I+1,J);
val(:) = rdx2;
tuples = [tuples; idp(:), idn(:), val(:) ];

idp = ind(I,J);
idn = ind(I,J-1);
val(:) = rdy2;
tuples = [tuples; idp(:), idn(:), val(:) ];

idp = ind(I,J);
idn = ind(I,J+1);
val(:) = rdy2;
tuples = [tuples; idp(:), idn(:), val(:) ];

idp = ind(I,J);
idn = ind(I,J);
val(:) = -2*rdx2 - 2*rdy2;
tuples = [tuples; idp(:), idn(:), val(:) ];


%
% treat boundary part
%
for j = 1:ny
for i = 1:nx
if i==1 | i==nx | j==1 | j==ny
    idx = ind(i,j);
    aa = 0;
    
    if i > 1
        cc = rdx2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(i-1,j), cc];
    elseif bctype(1,1) == bc_per
        cc = rdx2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(nx,j), cc];
        rb(idx) = rb(idx) - (bcval(1,2)-bcval(1,1))*cc;
    elseif bctype(1,1) == bc_neu
        rb(idx) = rb(idx) - bcval(1,1)*rdx;
    elseif bctype(1,1) == bc_dir
        cc = 2*rdx2;
        aa = aa - cc;
        rb(idx) = rb(idx) + bcval(1,1)*cc;
    end
    
    if i < nx
        cc = rdx2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(i+1,j), cc];
    elseif bctype(1,2) == bc_per
        cc = rdx2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(1,j), cc];
        rb(idx) = rb(idx) + (bcval(1,2)-bcval(1,1))*cc;
    elseif bctype(1,2) == bc_neu
        rb(idx) = rb(idx) + bcval(1,2)*rdx;
    elseif bctype(1,2) == bc_dir
        cc = 2*rdx2;
        aa = aa - cc;
        rb(idx) = rb(idx) + bcval(1,2)*cc;
    end
    
    if j > 1
        cc = rdy2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(i,j-1), cc ];
    elseif bctype(2,1) == bc_per
        cc = rdy2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(i,ny), cc ];
        rb(idx) = rb(idx) - (bcval(2,2)-bcval(2,1))*cc;
    elseif bctype(2,1) == bc_neu
        rb(idx) = rb(idx) - bcval(2,1)*rdy;
    elseif bctype(2,1) == bc_dir
        cc = 2*rdy2;
        aa = aa - cc;
        rb(idx) = rb(idx) + bcval(2,1)*cc;
    end
    
    if j < ny
        cc = rdy2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(i,j+1), cc ];
    elseif bctype(2,2) == bc_per
        cc = rdy2;
        aa = aa - cc;
        tuples(end+1,:) = [ idx, ind(i,1), cc ];
        rb(idx) = rb(idx) + (bcval(2,2)-bcval(2,1))*cc;
    elseif bctype(2,2) == bc_neu
        rb(idx) = rb(idx) + bcval(2,2)*rdy;
    elseif bctype(2,2) == bc_dir
        cc = 2*rdy2;
        aa = aa - cc;
        rb(idx) = rb(idx) + bcval(2,2)*cc;
    end
    
    % diag
    tuples(end+1,:) = [ idx, idx, aa ];
end
end
end

Lap = sparse(tuples(:,1), tuples(:,2), tuples(:,3), np,np,nbuf);

return
end


