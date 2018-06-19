function [v] = voigtEncode(m, off_scale)

if off_scale
	s = 2.0;
else
	s = 1.0;
end

v = [ m(1,1); m(2,2); s*m(1,2); m(3,3) ];


return
end

