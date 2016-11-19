

% analytical flux from internal boundary
% = (particle circumference) * (flux density)
flux_part_ana = sum(2*pi*partrad.*partdu)
if (1)
    % check numerical flux through external wall 
    % they should be equal
    flux_wall = 0;
    if bctype(2,1) == 1
        flux_wall = flux_wall - sum(bcval(2,1)-sol(:,1))/(dy/2)*dx;
    end
    if bctype(2,2) == 1
        flux_wall = flux_wall - sum(bcval(2,2)-sol(:,ny))/(dy/2)*dx;
    end
    flux_wall
end
if (1)
    % check numerical flux through internal particles
    % they should be equal to the numerical flux through external wall
    flux_part = 0;
    for j = 1:ny
    for i = 1:nx
    if tag(i,j) == 0
        if tag(i+1,j) == 1
            ff = (sol(i+1,j)-sol(i,j)) / dx * dy;
            flux_part = flux_part + ff;
        end
        if tag(i-1,j) == 1
            ff = (sol(i-1,j)-sol(i,j)) / dx * dy;
            flux_part = flux_part + ff;
        end
        if tag(i,j+1) == 1
            ff = (sol(i,j+1)-sol(i,j)) / dy * dx;
            flux_part = flux_part + ff;
        end
        if tag(i,j-1) == 1
            ff = (sol(i,j-1)-sol(i,j)) / dy * dx;
            flux_part = flux_part + ff;
        end
    end
    end
    end
    flux_part
end
if (1)
    % check source term
    dens_fluid = 0;
    sh = sinh(sol(:));
    sh(tag_nonfd) = 0;
    dens_fluid = sum(sh) * dx*dy * kappa2
end


