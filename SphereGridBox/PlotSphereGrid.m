function [fig] = PlotSphereFibonacciGrid ( ng, xg )

fig = figure ( );
clf
hold on

% draw sphere
[ x, y, z ] = sphere ( 40 );
shrink = 0.95;
x = shrink * x;
y = shrink * y;
z = shrink * z;
c = ones ( size ( z ) );
% surf ( x, y, z, c, 'EdgeColor','none');
surf ( x, y, z, c);

% draw points
plot3 ( xg(:,1), xg(:,2), xg(:,3), 'b.', 'Markersize', 20 );

%
if (1)
    for i = 1:ng
        if mod(i-1,2)==0
        text(xg(i,1),xg(i,2),xg(i,3), int2str(i), 'FontSize',12);
        end
    end
end

axis equal
grid on
view ( 3 )
xlabel ( '<--X-->' )
ylabel ( '<--Y-->' )
zlabel ( '<--Z-->' )
title ( sprintf ( '%d point Grid on Sphere', ng ), 'FontSize', 24 );

hold off

% print ( '-dpng', filename );
% print ( '-dpdf', filename );
% fprintf ( 1, '\n' );
% fprintf ( 1, '  Plot file saved to "%s".\n', filename );

return
end