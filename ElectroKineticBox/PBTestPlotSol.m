
%
% call FIGURE by yourself...
%

solsave = sol;
sol(tag_solid) = nan;

% imagesc(xcell,ycell,sol);
contourf(xcell,ycell,sol,-1.4:0.1:-0.2);
colorbar;
% contourf(xcell,ycell,phi);
% contourf(xcell,ycell,dphi);
% imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],owner');
% imagesc([xlo+dx/2,xhi-dx/2],[ylo+dy/2,yhi-dy/2],tag');
axis equal;
hold on;
contour(xcell,ycell,sdf,[0,0]);
hold off;

sol = solsave;
clear solsave;




