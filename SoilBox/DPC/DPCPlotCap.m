function DPCPlotCap(varargin)

DPCGlobals;

% figure(hfig);

% L
ll = el;
% R
rr = cr0;

% calculate X
xx = rr * (aa - ac.*exp(-ab.*ll)) + ll;

tt = linspace(0,pi/2,31);

xplot = ll + (xx-ll) * cos(tt);
yplot = (xx-ll)/rr * sin(tt);

plot(xplot,yplot,varargin{:});



return
end


