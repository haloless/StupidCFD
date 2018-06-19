
clear;

re = 2.0;

xs = linspace(0.0, 6.0, 7);
ys = linspace(0.0, 6.0, 7);
[xs,ys] = ndgrid(xs,ys);
xs = xs(:);
ys = ys(:);

figure;
plot(xs,ys,'o');
axis equal;


if 1
    xtest = [ 0.0, 0.0 ];
    % xtest = [ 1.0, 1.0 ];
    % xtest = [ 2.0, 2.0 ];
    % xtest = [ 3.1, 3.1 ];
    
    neigh = find(sqrt((xtest(1)-xs).^2 + (xtest(2)-ys).^2) < re);
    
    hold on;
    rectangle('Position',[xtest-re [re*2,re*2]], 'Curvature',[1,1]);
    plot(xtest(1),xtest(2),'xr', xs(neigh),ys(neigh),'+');
    hold off;
    
    [N,Nx,Ny] = mlsRegularShape2D(xtest, xs(neigh),ys(neigh), re);
    
    [sum(N), sum(xs(neigh).*N(:)), sum(ys(neigh).*N(:))]
    
    [Nx; Ny] * [xs(neigh), ys(neigh)]
    
end


