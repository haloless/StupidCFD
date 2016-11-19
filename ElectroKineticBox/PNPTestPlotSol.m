
%
% call FIGURE by yourself...
%

solsave = sol;

sol = reshape(sol,nx,ny);
sol(tag_solid) = nan;

if isempty(solrange)
    contourf(xcell,ycell,sol,'ShowText','on');
else
    contourf(xcell,ycell,sol,solrange, 'ShowText','on');
end
% imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],owner');
% imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],tag');
colorbar;
axis equal;
axis([xlo,xhi,ylo,3.5*D]);
% axis([7.5,12.5,7.5,12.5]);

hold on;
contour(xcell,ycell,sdf,[0,0],'LineColor','r');
hold off;

sol = solsave;
clear solsave;




