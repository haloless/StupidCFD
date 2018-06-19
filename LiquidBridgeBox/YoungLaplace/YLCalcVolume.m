function [vol] = YLCalcVolume(xs,rs)




rm = 0.5 * (rs(1:end-1) + rs(2:end));
xl = xs(2:end) - xs(1:end-1);
vol = sum(pi .* rm.^2 .* xl);



return
end

