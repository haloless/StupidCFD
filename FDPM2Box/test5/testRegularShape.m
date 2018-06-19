
clear;

re = 2.0;

% xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
% iprobe = 4;

xs = [0.0, 1.5, 3.0, 3.5, 4.0, 5.0, 6.0];
iprobe = 3;

ys = zeros(size(xs));
ys(iprobe) = 1;

% xx = xs(iprobe)-re:0.01:xs(iprobe)+re;
xx = 0:0.02:6.0;
ww = [];
for xint = xx
%for xint = [2.8]
    w = mlsRegularShape1D(xint, xs, re);
    ww(end+1) = w(iprobe);
    % ww(end+1,:) = w;
    % plot(xx,ww,'x-');
    % pause;
end

% figure;
plot(xx, ww,'-', xs,ys,'.r');
axis([0.0 6.0 -0.2 1.5]);


if 1
    % xtest = 2.;
    xtest = 3;
    [N,Nx] = mlsRegularShape1D(xtest, xs,re);
    
    dx = 1.0e-5;
    Nr = mlsRegularShape1D(xtest+dx, xs,re);
    Nl = mlsRegularShape1D(xtest-dx, xs,re);
    dN = (Nr-Nl) / (dx*2);
    
    Nx
    dN
    norm(Nx-dN)
    
    % Nx = dN;
    [sum(N)-1,  sum(xs.*N)-xtest, sum(Nx)-0, sum(xs.*Nx)-1]
    
end


