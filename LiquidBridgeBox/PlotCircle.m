function [] = PlotCircle(cen,rad,varargin)

pos = [cen(1)-rad,cen(2)-rad,rad*2,rad*2];
curv = [1,1];
rectangle('Position',pos,'Curvature',curv, varargin{:});


return
end

