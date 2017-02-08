function [] = PlotSDF(sdf, varargin)

contour(sdf.xg,sdf.yg,sdf.phig, varargin{:});

return
end

