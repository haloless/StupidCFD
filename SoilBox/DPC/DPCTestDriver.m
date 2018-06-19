
clear;

DPCGlobals;

sigx = 0; sigy = 0; sigz = 0; sigxz = 0; sigyz = 0; sigxy = 0;
epsx = 0; epsy = 0; epsz = 0; epsxz = 0; epsyz = 0; epsxy = 0;
depsx = 0; depsy = 0; depsz = 0; depsxz = 0; depsyz = 0; depsxy = 0;
mtype = 0;
ltype = 0;
tcut = 0; fcut = 0;
geop = 0;
nocon = 0;
istat = 0;
el = 0;

aki = 670; ak1 = 0.77; ak2 = 0.1;
agi = 50; ag1 = 0.77; ag2 = 0.433;
aa = 15; ab = 0.0025; ac = 14.9;
cr0 = 2; cr1 = 0; cr2 = 0;
aw = 0.15; ad = 0.5;
geop = 0;

global cons
global qgamma qdelta
global vv vh iux

sa = 6; 

vv = 8.0e-4;
qalpha = 3;
cons = 0.5;

nsteps = 4000;
iskip = 2;
iux = 0;
% layer = 0; % to have compacting material
layer = 1; % to have dilative material
istat = 1;

qgamma = 1/3 - 1/qalpha;
qdelta = 2/3 + 1/qalpha;

% compute initial cap
if 1
    DPCInitCap;
else
    % astart = 0;
    % astart = 100000;
    astart = 2;
end
disp(['astart=',num2str(astart)]);
% return


depsx = vv;
depsy = vv;
depsz = vv;
depsxz = 0;
depsyz = 0;
depsxy = 0;

p = geop;
s = zeros(5,1);

sigx = geop;
sigy = geop;
sigz = geop;
sigxz = 0;
sigyz = 0;
sigxy = 0;

el = astart;


if 1
    hfig = figure;
    ii = -2:0.05:3;
    ff = aa - ac.*exp(-ab.*ii);
    plot(ii,ff);
    % axis equal;
    % axis([-0.2,3,0,0.25]);
    axis([-2,3,0,0.25]);
    
    hold on;
    % ll = el;
    % rr = cr0;
    % xx = rr * (aa - ac.*exp(-ab.*ll)) + ll;
    % rectangle('Position',[ll*2-xx,0-(xx-ll)/rr,(xx-ll)*2,(xx-ll)/rr*2], 'Curvature',[1,1]);
    DPCPlotCap('-r');
    hold off;
    % return
end


disp('Begin isotropic consolidate');
data = [];
% load
ok = 0;
dx = vv;
x = 0;
for step = 1:1000
    xsave = x;
    psave = p;
    esave = el;
    x = x + dx;
    
    [p,s] = DPCModel(p,s);
    
    if p<=cons && abs(p-cons)<=8.0e-4
        % reached consolidation desired
        ok = 1;
    elseif p < cons
        ok = 0;
    else
        % overshoot, reduce increment
        dx = dx / 2;
        depsx = dx; 
        depsy = dx;
        depsz = dx;
        x = xsave;
        p = psave;
        el = esave;
        continue;
    end
    
    data(end+1,:) = [x,p];
    
    xpercent = x*100;
    disp(num2str([sigx,sigy,sigz,xpercent,xpercent,xpercent]));
    
    if 1
        figure(hfig);
        hold on;
        
        i1 = p * 3;
        sqrtj2 = sqrt(s(1)^2 + s(2)^2 + s(1)*s(2) + s(3)^2 + s(4)^2 + s(5)^2);
        plot(i1,sqrtj2,'x');
        
        % ll = el;
        % rr = cr0;
        % xx = rr * (aa - ac.*exp(-ab.*ll)) + ll;
        % rectangle('Position',[ll*2-xx,0-(xx-ll)/rr,(xx-ll)*2,(xx-ll)/rr*2], 'Curvature',[1,1]);
        DPCPlotCap('-k');
        
        hold off;
        
        pause;
    end
    
    if ok
        break;
    end
end
disp('End isotropic consolidate');

% after isotropic consolidation
ez0 = x;
epsx = x;
epsy = x;
epsz = x;


if 0
    % unload
    disp('Begin unload');
    ok = 0;
    cons = 0;
    % dx = -vv * 0.01;
    dx = -1e-5;
    depsx = dx;
    depsy = dx;
    depsz = dx;
    for step = 1:1000
        xsave = x;
        psave = p;
        esave = el;
        x = x + dx;
        
        [p,s] = DPCModel(p,s);
        
        if p>=cons && abs(p-cons)<=8.0e-4
            % reached consolidation desired
            ok = 1;
        elseif p > cons
            ok = 0;
        else
            % overshoot, reduce increment
            dx = dx / 2;
            depsx = dx; 
            depsy = dx;
            depsz = dx;
            x = xsave;
            p = psave;
            el = esave;
            continue;
        end
        
        data(end+1,:) = [x,p];
        
        xpercent = x*100;
        disp(num2str([sigx,sigy,sigz,xpercent,xpercent,xpercent]));
        
        if ok
            break;
        end
    end
    disp('End unload');
end



if 1
    figure;
    plot(data(:,1)*3*100,data(:,2),'x-');
    xlabel('volumetric strain (%)');
    ylabel('mean normal stress');
end

if 1
    disp('Begin uniaxial test');
    
    iprint = 1;
    kount = 0;
    ikcount = 0;
    
    vv = 1.0e-4;
    
    ccgk = 1.5 * aki*(1-ak1) / (agi*(1-ag1));
    anu = (ccgk-1) / (2*ccgk+1);
    vh = -anu * vv;
    
    if 0
        % strain controled
        iux = 1;
        vh = 0;
    else
        % stress controled
        iux = 0;
    end
    
    nunld = 0;
    
    for kount = 1:nsteps
        
        [p,s] = DPCLoadTriaxial(nunld,p,s);
        
        
        epsx = epsx + depsx;
        epsy = epsy + depsy;
        epsz = epsz + depsz;
        
        vepsz = (epsz-ez0) * 100;
        
        if mod(kount,2) == 0
            disp(num2str([sigx,sigy,sigz,epsx*100,epsy*100,epsz*100]));
        end
        
        if mod(kount,2) == 0
            figure(hfig);
            plot(ii,ff);
            % axis([-0.2,3,0,0.25]);
            axis([-2,3,0,0.25]);

            
            hold on;
            
            i1 = p * 3;
            sqrtj2 = sqrt(s(1)^2 + s(2)^2 + s(1)*s(2) + s(3)^2 + s(4)^2 + s(5)^2);
            plot(i1,sqrtj2,'x');
            
            % ll = el;
            % rr = cr0;
            % xx = rr * (aa - ac.*exp(-ab.*ll)) + ll;
            % rectangle('Position',[ll*2-xx,0-(xx-ll)/rr,(xx-ll)*2,(xx-ll)/rr*2], 'Curvature',[1,1]);
            DPCPlotCap('-k');
            
            hold off;
            
            pause;
        end
        
        if vepsz>sa
            break;
        end
    end
    
    disp('End uniaxial test');

end


