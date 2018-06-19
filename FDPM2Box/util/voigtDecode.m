function [m] = voigtDecode(v, off_scale)

if off_scale
	s = 0.5;
else
	s = 1.0;
end

m = [ v(1), s*v(3), 0; s*v(3), v(2), 0; 0, 0, v(4) ];


return
end

