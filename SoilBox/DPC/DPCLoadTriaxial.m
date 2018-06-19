function [p,s] = DPCLoadTriaxial(nunld,p,s)

DPCGlobals;

global cons
global qgamma qdelta
global vv vh iux

err = 1e-3;

pa = p;
sa = s(1);
sb = s(2);
sc = s(3);
sd = s(4);
se = s(5);
esave = el;
btrial = vh;

% set strain incr
strain(nunld);

% solve model
[p,s] = DPCModel(p,s);


if iux == 1
    % strain controlled model
    % this is done
else
    % stress controlled model
    % adjust strain to achieve stress state
    
    sigxx = p + s(1);
    sigzz = p - s(1) - s(2);
    test = sigzz*qgamma + sigxx*qdelta - cons;
    if abs(test) < err
        return;
    end
    
    tsave = test;
    if test < 0
        if btrial < 0
            db = btrial;
        elseif btrial >= 0
            db = -btrial/2;
        end
    elseif test > 0
        if btrial < 0
            db = -btrial/2;
        elseif btrial >= 0
            db = btrial;
        end
    end
    
    % 60
    bt = btrial + db;
    
    ok = 0;
    niter = 500;
    for iter = 1:niter
        vh = bt;
        p = pa;
        s = [ sa, sb, sc, sd, se ]';
        el = esave;
        
        strain(nunld);
        [p,s] = DPCModel(p,s);
        
        sigxx = p + s(1);
        sigzz = p - s(1) - s(2);
        test = sigzz*qgamma + sigxx*qdelta - cons;
        if abs(test) < err
            ok = 1; break;
        end
        
        btrial = bt;
        bt = bt - test*db/(test-tsave);
        db = bt - btrial;
        tsave = test;
    end
    
    if ~ok
        error('stress control failed');
    end
end



return
end


function strain(nunld)
DPCGlobals;
global vv vh iux

depsx = vh;
depsy = vh;
depsz = vv;
depsxz = 0;
depsyz = 0;
depsxy = 0;

if nunld
    depsx = -depsx;
    depsy = -depsy;
    depsz = -depsz;
end

return
end

