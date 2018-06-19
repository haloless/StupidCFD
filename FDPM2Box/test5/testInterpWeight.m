
clear;

xs = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
% xs = [0.0, 1.5, 3.0, 3.5, 4.0, 5.0, 6.0];

% re = 2.0;
re = 2.5;
% re = 3.0;

xx = 1.0:0.01:5.0;
ww = [];
for xint = xx
    w = mlsInterpWeight1D(xint, xs, re);
    ww(end+1) = w(4);
    % ww(end+1,:) = w;
    % plot(xx,ww,'x-');
    % pause;
end


% plot(ww,'-');
plot(xx,ww,'-'); %axis([1 5 0 1])


if 1
    % use finite difference to check weight derivatives
    
    dx = 1.0e-5;
    
    xtest = 4.1;
    [w,wx] = mlsInterpWeight1D(xtest, xs, re);
    
    wp = mlsInterpWeight1D(xtest+dx, xs, re);
    wm = mlsInterpWeight1D(xtest-dx, xs, re);
    dw = (wp-wm) ./ (dx*2);
    
    wx
    dw
    norm(wx-dw)
end


