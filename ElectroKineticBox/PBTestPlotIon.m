
%
% call FIGURE by yourself...
%

zz = 1.0;
sol_cation = exp(-zz .* sol);
sol_anion = exp(-(-zz) .* sol);
sol_cation(tag_solid) = nan;
sol_anion(tag_solid) = nan;


figure;
contourf(xcell,ycell,sol_cation);
colorbar;
title('cation');
axis equal;
hold on;
contour(xcell,ycell,sdf,[0,0]);
hold off;

figure;
contourf(xcell,ycell,sol_anion);
colorbar;
title('anion');
axis equal;
hold on;
contour(xcell,ycell,sdf,[0,0]);
hold off;




