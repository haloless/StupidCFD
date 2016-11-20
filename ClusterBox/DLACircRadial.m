
function [] = DLACircRadial(np)

% np = 32000;

xlo = -100;
xhi = 100;
ylo = -100;
yhi = 100;

a = 1.0;
d = a * 2;

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
        spd = 0.1;
        ang = rand(1) * 2*pi;
        
        while 1
            ang = ang + rand(1)-0.5;
            
            u = cos(ang) * spd;
            v = sin(ang) * spd;
            
            x = x + u;
            y = y + v;
            
            if x^2 + y^2 > rout^2
                break;
            end
            
            if check_col(a, x,y, cnt,xs,ys) == 1
                ok = 1;
                break;
            end
        end
    end
    
	% while check_col(x,y,nx,ny,mat) == 0
		% [x,y,spd,ang] = move_part(x,y,spd,ang,nx,ny);
	% end
	
	% i = round(x);
	% j = round(y);
	% if i>=1 && i<=nx && j>=1 && j<=ny
		% mat(i,j) = 1;
	% end
	
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
prompt = ['part=',int2str(cnt), ';Rg=',num2str(rg)];
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


function [col] = check_col(a, x,y, n,xs,ys)
    d2 = 4*a^2;
    
    col = 0;
    for i = 1:n
        dist2 = (x-xs(i))^2 + (y-ys(i))^2;
        if dist2 < d2
            col = 1;
            break;
        end
    end
    
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


