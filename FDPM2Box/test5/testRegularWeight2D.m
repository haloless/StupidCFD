
clear;

re = 2.0;
% re = 2.5;
% re = 3.0;


xs = linspace(0.0, 6.0, 7);
ys = linspace(0.0, 6.0, 7);


% xs = [0.0, 1.5, 3.0, 3.5, 4.0, 5.0, 6.0];


if 1
    % use finite difference to check weight derivatives
    
    dx = 1.0e-5;
    
    xtest = [2.4, 2.0];
    [w,wx,wy] = mlsRegularWeight2D(xtest, xs,ys, re);
    
    wp = mlsRegularWeight2D(xtest+[dx,0], xs,ys, re);
    wm = mlsRegularWeight2D(xtest-[dx,0], xs,ys, re);
    dwx = (wp-wm) ./ (dx*2);
    [wx;dwx]
    
    wp = mlsRegularWeight2D(xtest+[0,dx], xs,ys, re);
    wm = mlsRegularWeight2D(xtest-[0,dx], xs,ys, re);
    dwy = (wp-wm) ./ (dx*2);
    [wy;dwy]
end


