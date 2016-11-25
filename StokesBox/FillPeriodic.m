
function [uper] = FillPeriodic(u,nx,ny)

uper = zeros(nx+1,ny+1);
uper(1:nx,1:ny) = u;
uper(:,ny+1) = uper(:,1);
uper(nx+1,:) = uper(1,:);

return
end
