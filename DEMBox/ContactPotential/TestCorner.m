
clear all;

corner = MakeCorner([0;0],[0;3],[3;0]);
% corner = MakeCorner([0;0],[-0.5;1],[1;0.5]);
% corner = MakeCorner([0;0],[-0.5;-1],[-0.5;1]);

hfig = figure;
hold on;
plot([corner.xa(1),corner.xw(1),corner.xb(1)],[corner.xa(2),corner.xw(2),corner.xb(2)],'kx-');
if 0
    for x = -1:0.1:1
    for y = -1:0.1:1
        phi = CornerPotential(corner,x,y);
        if phi >= 0
            plot(x,y,'xb');
        else
            plot(x,y,'or');
        end
    end
    end
end
hold off;
axis equal;
axis([-0.5 3 -0.5 3]);

xbase = [0.6;0.5];
shape = MakeSuperEllipse(1.0,0.5, 2,2, xbase(1),xbase(2), 0/180*pi);

dh = 0.25;
% xbase = [0.5; 0.5];
% xbase = [0.5; 0.3];
% xbase = [0.5; 0.1];
% xbase = [0.58; 0.1];

dw = CornerPotential(corner,xbase(1),xbase(2));
[gw,hw] = CornerDeriv(corner,xbase(1),xbase(2),dh);

wfunc = @(x,y) CornerApprox(x,y, dw,gw,hw,xbase);

figure(hfig);
% figure;
hold on;
PlotShape(shape,{'m'});
ezplot(wfunc);
plot(xbase(1),xbase(2),'rx');
hold off;

