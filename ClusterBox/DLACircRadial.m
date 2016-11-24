
function [] = DLACircRadial(np,stickprob)

% np = 32000;

% np = 100;

twopi = pi * 2;

xlo = -100;
xhi = 100;
ylo = -100;
yhi = 100;

a = 1.0;
d = a * 2;

if ~exist('stickprob','var')
	stickprob = 1.0;
end

% stickprob = 1;
% stickprob = 0.1;
% stickprob = 0.01;


figure;

xs = zeros(np,1);
ys = zeros(np,1);

xs(1) = 0;
ys(1) = 0;
phi = rand(1) * 2*pi;
xs(2) = d * cos(phi);
ys(2) = d * sin(phi);

cnt = 2;


while cnt < np
    
    [x0,y0] = calc_cen(cnt,xs,ys);
    rg = calc_gyro(cnt,xs,ys);
    rin = rg * 4.0;
    rout = rg * 5.0;
    
    ok = 0;
    
    while ok == 0
        phi = rand(1) * 2*pi;
        x = rin * cos(phi);
        y = rin * sin(phi);
        spd = 1.0 * a;
        ang = rand(1) * 2*pi;
        
        while 1
            % ang = ang + rand(1)-0.5;
            ang = rand(1) * 2*pi;
            
			eu = cos(ang);
			ev = sin(ang);
            u = eu * spd;
            v = ev * spd;
            
            x = x + u;
            y = y + v;
            
            if x^2 + y^2 > rout^2
				% has gone too far, restart random walk
                break;
            end
            
			x = x - u;
			y = y - v;
			
			[isec,dsec] = intersect_circ(d, x,y,eu,ev,spd, cnt,xs,ys);
			if isec > 0
				x = x + eu * dsec;
				y = y + ev * dsec;
				icol = 1;
				dist = d;
				ex = (x-xs(icol)) / dist;
				ey = (y-ys(icol)) / dist;
			else
				x = x + u;
				y = y + v;
				icol = 0;
			end
			
			% [icol,dist,ex,ey] = check_col(a, x,y, cnt,xs,ys);
			
            if icol > 0
				stick = rand(1);
				if stick < stickprob
					% can stick to the current site
					ok = 1;
					break;
				else
					% cannot stick, rebound
					% veln = u*ex + v*ey;
					% u2 = -veln*ex + (u-veln*ex);
					% v2 = -veln*ey + (v-veln*ey);
					% u2 = -veln*ex;
					% v2 = -veln*ey;
					% veln = (d - dist) * 10;
					% u2 = veln * ex;
					% v2 = veln * ey;
					u2 = -eu*dsec;
					v2 = -ev*dsec;
					x = x + u2;
					y = y + v2;
					% ang = rand(1)*twopi;
					% ang = ang + pi;
					% ang = atan2(ey,ex);
				end
			else
				% ang = ang + rand(1)-0.5;
				% ang = rand(1) * 2*pi;
			end
        end
    end
    
	if (1)
		
	end
	
	
    cnt = cnt + 1;
    xs(cnt) = x;
    ys(cnt) = y;
    
	if mod(cnt,100)==0 ||1
        prompt = ['part=',int2str(cnt)];
        plot(x0,y0,'xr');
		hold on;
        for i = 1:cnt
            rectangle('Position',[xs(i)-a,ys(i)-a,d,d], 'Curvature',[1,1]);
        end
        rectangle('Position',[x0-rg,y0-rg,rg*2,rg*2], 'Curvature',[1,1]);
        % rectangle('Position',[x0-rin,y0-rin,rin*2,rin*2], 'Curvature',[1,1]);
        % rectangle('Position',[x0-rout,y0-rout,rout*2,rout*2], 'Curvature',[1,1]);
        hold off;
        axis equal;
        axis([xlo xhi ylo yhi]);
        title(prompt);
        drawnow;
        % pause;
	end
end

rg = calc_gyro(cnt,xs,ys);

figure;
prompt = ['part=',int2str(cnt), ';stick=',num2str(stickprob), ';Rg=',num2str(rg)];
hold on;
for i = 1:cnt
    rectangle('Position',[xs(i)-a,ys(i)-a,d,d], 'Curvature',[1,1]);
end
hold off;
axis equal;
axis([xlo xhi ylo yhi]);
title(prompt);
drawnow;


return
end

function [x0,y0] = calc_cen(n,xs,ys)
    x0 = mean(xs(1:n));
    y0 = mean(ys(1:n));
return
end

function [rg] = calc_gyro(n,xs,ys)
    [x0,y0] = calc_cen(n,xs,ys);
    rg2 = 1.0/n * sum((xs(1:n)-x0).^2 + (ys(1:n)-y0).^2);
    rg = sqrt(rg2);
return
end


function [col,dist,ex,ey] = check_col(a, x,y, n,xs,ys)
    d2 = 4*a^2;
    
	col = 0;
	dist = 0;
	ex = 0;
	ey = 0;
	
	if 1
		
		I = 1:n;
		r2 = (x-xs(I)).^2 + (y-ys(I)).^2;
		[r2min,imin] = min(r2);
		
		if r2min <= d2
			col = imin;
			dist = sqrt(r2min);
			ex = (x-xs(imin)) / dist;
			ey = (y-ys(imin)) / dist;
		end
	end
	
    % for i = 1:n
        % dist2 = (x-xs(i))^2 + (y-ys(i))^2;
        % if dist2 <= d2
            % col = i;
			% dist = sqrt(dist2);
			% ex = (x-xs(i)) / dist;
			% ey = (y-ys(i)) / dist;
            % break;
        % end
    % end
    
	% [isec,dsec] = intersect_circ(a*2, x,y,
	
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


