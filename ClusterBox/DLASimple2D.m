
function [mat] = mydla(np)

nx = 256;
ny = 256;

mat = zeros(nx,ny);
mat(:,ny) = 1;

np = 32000;

figure;

for ip = 1:np
	
	[x,y,spd,ang] = gen_part(nx,ny);
	
	while check_col(x,y,nx,ny,mat) == 0
		[x,y,spd,ang] = move_part(x,y,spd,ang,nx,ny);
	end
	
	i = round(x);
	j = round(y);
	if i>=1 && i<=nx && j>=1 && j<=ny
		mat(i,j) = 1;
	end
	
	if mod(ip,100)==0 || ip==np
		prompt = ['part=',int2str(ip)];
		imshow(mat);
		title(prompt);
		drawnow;
	end
end

return
end

function [col] = is_col(i,j,nx,ny,mat)
	if i<1 || i>nx || j<1 || j>ny
		col = 0;
	else
		col = (mat(i,j)==1);
	end
return
end

function [col] = check_col(x,y,nx,ny,mat)
	i = round(x);
	j = round(y);
	col = ...
	is_col(i-1,j-1,nx,ny,mat) || is_col(i,j-1,nx,ny,mat) || is_col(i+1,j-1,nx,ny,mat) || ...
	is_col(i-1,j,nx,ny,mat) || is_col(i,j,nx,ny,mat) || is_col(i+1,j,nx,ny,mat) || ...
	is_col(i-1,j+1,nx,ny,mat) || is_col(i,j+1,nx,ny,mat) || is_col(i+1,j+1,nx,ny,mat);
	
return
end

function [x,y,spd,ang] = gen_part(nx,ny)
	x = rand(1) * nx;
	y = 0;
	spd = 1.0;
	ang = rand(1) * 2*pi;
return
end

function [x,y,spd,ang] = move_part(x,y,spd,ang,nx,ny)
	ang = ang + rand(1)-0.5;
	
	if ang > pi
		% ang = pi;
	elseif ang < 0
		% ang = 0;
	end
	
	u = cos(ang) * spd;
	v = sin(ang) * spd;
	if v < 0
		v = -v;
	end
	
	x = x + u;
	y = y + v;
	
	if x > nx
		x = x - nx;
	elseif x <= 0
		x = x + nx;
	end
	if y > ny
		y = y - ny;
	elseif y <=0 
		y = y + ny;
	end
return
end


