
clear all;

D = 1.0;
% Lx = 10.0;
% Ly = 10.0;
% Lx = 15.0;
% Ly = 15.0;
Lx = 20.0;
Ly = 20.0;
xlo = 0.0; xhi = Lx;
ylo = 0.0; yhi = Ly;
xcen = 0.5 * (xlo+xhi);
ycen = 0.5 * (ylo+yhi);

% nx = 128;
% ny = 128;
% dx = Lx / nx; 
% dy = Ly / ny;
dx = D / 12.8;
dy = D / 12.8;
% dx = D / 8;
% dy = D / 8;
% dx = D / 10;
% dy = D / 10;
nx = round(Lx/dx);
ny = round(Ly/dy);
dx2 = dx*dx;
dy2 = dy*dy;

% 
dl = 0.5 * min(dx,dy);


ncell = nx * ny;
xcell = linspace(xlo+dx*0.5,xhi-dx*0.5,nx);
ycell = linspace(ylo+dy*0.5,yhi-dy*0.5,ny);
[xcell,ycell] = ndgrid(xcell,ycell);

% conductivity
kf = 1.0;
ks = 0.0;
kr = ks / kf;

% BC
uylo = 0.0;
uyhi = 1.0 * Ly;
ucen = 0.0;

%
sdf = zeros(nx,ny);
nvecx = zeros(nx,ny);
nvecy = zeros(nx,ny);
for j = 1:ny
for i = 1:nx
    xx = xcell(i,j) - xcen;
    yy = ycell(i,j) - ycen;
    rr = sqrt(xx^2 + yy^2);
    sdf(i,j) = rr - D/2;
    nvecx(i,j) = xx / max(rr,D*1.0e-4);
    nvecy(i,j) = yy / max(rr,D*1.0e-4);
end
end

%
tag = SharpIBTag(nx,ny,sdf);
isolid = find(tag<=0);
ifluid = find(tag>0);
inull = find(tag==-1);
ifixed = find(tag==0);

% figure;
% imagesc(tag);

%
bctype = [ 2, 2; 1, 1 ];
bcvalue = [ 0, 0; uylo, uyhi ];
[Lap,rb] = LaplacianOp2D(nx,ny,dx,dy, bctype,bcvalue);
disp(['Created LapOp']);

%
MInt = sparse(ncell,ncell);
MInt1 = sparse(ncell,ncell);
MInt2 = sparse(ncell,ncell);
for j = 1:ny
for i = 1:nx
    if tag(i,j) == 0
        ind = i + (j-1)*nx;
        
        dist = abs(sdf(i,j));
        distx = dist * nvecx(i,j);
        disty = dist * nvecy(i,j);
        
        % ghost point position
        xx = xcell(i,j); 
        yy = ycell(i,j);
        
        % surface position
        xb = xx + distx;
        yb = yy + disty;
        
        % mirror point
        xm = xb + distx;
        ym = yb + disty;
        
        [ ws, is, js ] = GridInterpCoef2D(xm,ym, xlo,ylo,dx,dy);
        for k = 1:4
            indk = is(k) + (js(k)-1)*nx;
            MInt(ind,indk) = ws(k);
            
            if tag(is(k),js(k)) == -1
                error('solid point included');
            end
        end
        
        
        % extended point 1
        x1 = xb + dl * nvecx(i,j);
        y1 = yb + dl * nvecy(i,j);
        
        [ ws, is, js ] = GridInterpCoef2D(x1,y1, xlo,ylo,dx,dy);
        for k = 1:4
            indk = is(k) + (js(k)-1)*nx;
            MInt1(ind,indk) = ws(k);
            
            if tag(is(k),js(k)) == -1
                error('solid point included');
            end
        end
        
        % extended point 2
        x2 = xb + 2*dl * nvecx(i,j);
        y2 = yb + 2*dl * nvecy(i,j);
        
        [ ws, is, js ] = GridInterpCoef2D(x2,y2, xlo,ylo,dx,dy);
        for k = 1:4
            indk = is(k) + (js(k)-1)*nx;
            MInt2(ind,indk) = ws(k);
            
            if tag(is(k),js(k)) == -1
                error('solid point included');
            end
        end
    end
end
end
disp(['Created IntOp']);


Dfixed = zeros(ncell,1);
Dfixed(ifixed) = 1;
Dfixed = spdiags(Dfixed,0,ncell,ncell);

uvec = zeros(ncell,1);
for iter = 1:1000
    
    %
    uint = uvec;
    
    if (0)
        % Dirichlet BC
        for fixed_iter = 1:200
            usave = uint;
            utmp = MInt * uint;
            uint(ifixed) = 2*ucen - utmp(ifixed);
            
            rinc = norm(uint - usave);
            if fixed_iter == 1
                rfixed0 = rinc;
            elseif rinc<rfixed0*1e-8 | rinc<1.0e-8
                disp(['fixed-point converged: iter=',int2str(fixed_iter)]);
                break;
            end
            % disp(['|uinc|=', num2str(rinc)]);
        end
    end
    if (1)
        % Neumann BC
        for fixed_iter = 1:1
            usave = uint;
            
            utmp1 = MInt1 * uint;
            utmp2 = MInt2 * uint;
            utmp3 = -1.0/3.0 * (2*dl*ucen - 4.0 .* utmp1(ifixed) + utmp2(ifixed));
            
            utmp = MInt * uint;
            uint(ifixed) = 2*utmp3 - utmp(ifixed);
            
            rinc = norm(uint - usave);
            if fixed_iter == 1
                rfixed0 = rinc;
            elseif rinc<rfixed0*1e-8 | rinc<1.0e-8
                disp(['fixed-point converged: iter=',int2str(fixed_iter)]);
                break;
            end
            % disp(['|uinc|=', num2str(rinc)]);
        end
    end
    if (0)
        % Neumann BC2
        udist = abs(sdf(ifixed));
        
        for fixed_iter = 1:200
            usave = uint;
            
            utmp = MInt * uint;
            uint(ifixed) = utmp(ifixed) - 2*ucen.*udist;
            
            rinc = norm(uint - usave);
            if fixed_iter == 1
                rfixed0 = rinc;
            elseif rinc<rfixed0*1e-8 | rinc<1.0e-8
                disp(['fixed-point converged: iter=',int2str(fixed_iter)]);
                break;
            end
            % disp(['|uinc|=', num2str(rinc)]);
        end
    end
    % break;
    
    ularge = 1.0e12;
    
    A = Lap + ularge.*Dfixed;
    b = -rb + ularge.*(Dfixed*uint);
    x = A \ b;
    
    uprev = uvec;
    uvec = x;
    
    udiff = uvec - uprev;
    udiff(isolid) = 0;
    rdiff = norm(udiff);
    
    disp(['iter=',int2str(iter), '; |udiff|=',num2str(rdiff)]);
    if rdiff < 1.0e-6
        disp(['converged']);
        break;
    end
end


% sol = uvec;
sol = uint;
sol = reshape(sol, nx,ny);
% sol(inull) = NaN;
figure;
contourf(xcell,ycell,sol,32);
axis equal;
hold on;
contour(xcell,ycell,sdf,[0,0]);
hold off;

fyhi = sum((uyhi - sol(:,ny)) / (dy*0.5)) / nx
fylo = -sum((uylo - sol(:,1)) / (dy*0.5)) / nx



