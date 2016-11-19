
function [ Aout,rout ] = ModifyLap2D(Ain,rin,nx,ny,dx,dy,xlo,ylo,sdf,tag,owner,npart,partbc,partval)
% L(x) = A.x + r

bc_neu = 0;
bc_dir = 1;

Aout = Ain;
rout = rin;

% control point position
dh = max(dx,dy);
len = dh * 1.5;

ind = reshape(1:nx*ny, nx,ny);

for j = 1:ny
for i = 1:nx
    % ghost cell need modification
    if tag(i,j) == 0
        
        idx = ind(i,j);
        
        % distance and normal from wall
        dist = abs(sdf(i,j));
        nvecx = (sdf(i+1,j)-sdf(i-1,j)) / (dx*2);
        nvecy = (sdf(i,j+1)-sdf(i,j-1)) / (dy*2);
        nnorm = sqrt(nvecx^2 + nvecy^2);
        nvecx = nvecx / nnorm;
        nvecy = nvecy / nnorm;
        
        % this point
        % xx = xcell(i,j);
        % yy = ycell(i,j);
        xx = xlo + (i-0.5)*dx;
        yy = ylo + (j-0.5)*dy;
        
        % on boundary
        xb = xx + dist*nvecx;
        yb = yy + dist*nvecy;
        
        % 1st control point
        l1 = len * 1;
        x1 = xb + l1*nvecx;
        y1 = yb + l1*nvecy;
        % 2nd control point
        l2 = len * 2;
        x2 = xb + l2*nvecx;
        y2 = yb + l2*nvecy;
        
        % which particle
        ipart = owner(i,j);
        bctype = partbc(ipart);
        bcval = partval(ipart);
        
        if bctype == bc_neu
            % Neumann type, use 2 control points
            c0 = -(l1+l2) / (dist+l1) / (dist+l2);
            c1 = (dist-l2) / (dist+l1) / (l1-l2);
            c2 = (dist-l1) / (dist+l2) / (l2-l1);
            
            Aout(idx,:) = 0;
            Aout(idx,idx) = c0;
            rout(idx) = -bcval;
            
            [ is, js, ws ] = CellInterpCoef2D(x1,y1, nx,ny,dx,dy,xlo,ylo);
            for k = 1:4
                idk = ind(is(k),js(k));
                if tag(idk) ~= 1
                    error('bad interp');
                end
                Aout(idx,idk) = Aout(idx,idk) + c1*ws(k);
            end
            [ is, js, ws ] = CellInterpCoef2D(x2,y2, nx,ny,dx,dy,xlo,ylo);
            for k = 1:4
                idk = ind(is(k),js(k));
                if tag(idk) ~= 1
                    error('bad interp');
                end
                Aout(idx,idk) = Aout(idx,idk) + c2*ws(k);
            end
        elseif bctype==bc_dir
            % Dirichlet type, use 1 control point
            c0 = l1 / (l1+dist);
            c1 = dist / (l1+dist);
            
            Aout(idx,:) = 0;
            Aout(idx,idx) = c0;
            rout(idx) = -bcval;
            
            [ is, js, ws ] = CellInterpCoef2D(x1,y1, nx,ny,dx,dy,xlo,ylo);
            for k = 1:4
                idk = ind(is(k),js(k));
                if tag(idk) ~= 1
                    error('bad interp');
                end
                Aout(idx,idk) = Aout(idx,idk) + c1*ws(k);
            end
        end
    end
end
end



return
end

