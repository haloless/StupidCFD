function DPCMatLaw()

DPCGlobals;

conv = 1.0e-3;
niter = 60;


% volumetric and deviatoric
dev = depsx + depsy + depsz;
devo3 = dev / 3;
dexx = depsx - devo3;
deyy = depsy - devo3;

% initial stress
press = (sigx + sigy + sigz) / 3;
sxx = sigx - press;
syy = sigy - press;
sxz = sigxz;
syz = sigyz;
sxy = sigxy;
% initial invariant
sj1i = press * 3;
sj2i = sqrt(sxx^2 + syy^2 + sxx*syy + sxz^2 + syz^2 + sxy^2);
% initial
capl = BIGL(ep);
xl = X(ep);
evpi = EVP(xl);

% elastic
threek = 3 * BMOD(sj1i);
twog = 2 * SMOD(sj2i);

% elastic trial
sj1 = threek * dev + sj1i;
sxx = sxx + twog * dexx;
syy = syy + twog * deyy;
sxz = sxz + twog * depsxz;
syz = syz + twog * depsyz;
sxy = sxy + twog * depsxy;
ratio = 1.0;
mtype = 1;

%
tencut = min(fcut,tcut+3*geop);
if sj1 <= tencut
    % tensile coding
    sj1 = tencut;
    mtype = 0;
    
    if ltype~=2 && ep>astart
        % tension dilatency coding
        ell = max(astart, ep+conv*F(ep));
        xll = X(ell);
        denom = EVP(xll) - evpi;
        if denom > 0
            devp = dev - (sj1-sj1i)/threek;
            ep = ep + devp*(ell-ep)/denom;
            ep = min(astart,ep);
        else
            ep = astart;
        end
    end
    
else
    % 10 check failure envelop
    sj2 = sqrt(sxx^2 + syy^2 + sxx*syy + sxz^2 + syz^2 + sxy^2);
    if sj1 <= capl
        %
        tmises = SJ2C(capl,xl,capl);
        fj1 = F(sj1);
        ff = sj2 - min(fj1,tmises);
        if ff > 0
            % failure envelop violated
            % yield surface calculation
            mtype = 2;
            dfdj1 = 0;
            if fj1 < tmises
                dfdj1 = (F(sj1+conv*sj2)-fj1)/(conv*sj2);
            end
            dlamd = ff / (3*threek*dfdj1^2 + 0.5*twog);
            devp = -3 * dfdj1 * dlamd;
            sj1e = sj1;
            sj1 = sj1 - threek * devp;
            
            % now returned back to the failure envelop
            % dilatancy and corner coding
            if ltype==1 && ep>astart % compacting material
                % 60
                ell = min(sj1e, ep+conv*sj2);
                xll = X(ell);
                denom = EVP(xll) - evpi;
                if denom > 0
                    % 65
                    dilat = (ell-ep) / denom;
                    elint = ep;
                    ep = max(astart, elint+devp*dilat);
                    if sj1 > BIGL(ep)
                        sj1 = (dilat*sj1e + threek*elint) / (dilat+threek);
                        ep = max(astart, sj1);
                    end
                else % dilative material
                    sj1 = sj1e;
                    ep = max(astart, sj1);
                end
            else
                sj1 = min(sj1, capl);
            end
            
            % 70
            fj1 = F(sj1);
            ratio = min(fj1,tmises) / sj2;
            
        end
    else
        % cap hardening 
        % 30 
        sj2check = SJ2C(sj1,xl,capl);
        if sj1>xl || sj2>sj2check
            % 40
            sj1e = sj1;
            sj2e = sj2;
            
            ell = ep;
            elr = sj1e;
            if sj1e >= xl
                fl = (ep-sj1e) / (ep-xl);
            else
                fl = 2 * sj2e / (sj2e+SJ2C(sj1e,xl,capl)) - 1;
            end
            
            xr = X(elr);
            sj1r = sj1e - threek*(EVP(xr)-evpi);
            fr = (xr-sj1r) / (elr-xr);
            
            comp = conv * F((fl*xr-fr*xl)/(fl-fr));
            if abs(sj1)+sj2 >= comp
                mtype = 3;
                fold = 0;
                ok = 0;
                for iter = 1:niter
                    ep = (fl*elr-fr*ell) / (fl-fr);
                    xl = X(ep);
                    devp = EVP(xl) - evpi;
                    sj1 = sj1e - threek*devp;
                    capl = BIGL(ep);
                    if sj1 >= xl
                        fc = (ep-sj1) / (ep-xl);
                    end
                    if sj1 <= capl
                        fc = (xl-sj1) / (capl-xl);
                    end
                    
                    if sj1<xl && sj1>capl
                        sj2 = SJ2C(sj1,xl,capl);
                        delj1 = conv * (xl-sj1);
                        desp = 0;
                        if abs(delj1) > 1.0e-9
                            desp = (devp/6) * (delj1/(sj2-SJ2C(sj1+delj1,xl,capl)));
                        end
                        
                        sj2try = sj2 + twog*desp;
                        fc = (sj2e-sj2try) / (sj2e+sj2try);
                        
                        if abs(sj2e-sj2try) <= comp
                            ok = 1; break;
                        end
                        if fc>0 && (sj1-capl)<=delj1
                            ok = 1; break;
                        end
                    end
                    
                    % 300
                    if fc <= 0
                        elr = ep;
                        fr = fc;
                        if fold < 0
                            fl = 0.5 * fl;
                        end
                    else
                        ell = ep;
                        fl = fc;
                        if fold > 0
                            fr = 0.5 * fr;
                        end
                    end
                    
                    % 190
                    fold = fc;
                end
                
                if ~ok
                    warning('Hardening failed to converge');
                    nocon = 1;
                    sj1 = min(sj1,xl);
                    if sj1 < BIGL(elr)
                        sj1 = capl;
                    end
                    sj2 = min(sj2e, SJ2C(sj1,xl,capl));
                else
                    % disp(['Hardening converged, iter=',int2str(iter)]);
                end
                
                % 195
                ratio = 0;
                if abs(sj2e) > 1.0e-9
                    ratio = sj2 / sj2e;
                end
            end
        end
        
        
    end
end


% 200
sxx = sxx * ratio;
syy = syy * ratio;
sxz = sxz * ratio;
syz = syz * ratio;
sxy = sxy * ratio;
press = sj1 / 3;
sigx = sxx + press;
sigy = syy + press;
sigz = press - sxx - syy;
sigxz = sxz;
sigyz = syz;
sigxy = sxy;

return
end

function [ret] = F(i1) 
% envelop: A - C * exp(-B*I1) 
DPCGlobals;
ret = aa - ac*exp(-ab*i1);
end

function [ret] = BIGL(x)
% cap L = max(l,Linit)
DPCGlobals;
ret = max(astart,x);
end

function [ret] = R(capl)
% cap R(L)
DPCGlobals;
ret = cr0 * (1-cr1*exp(-cr2*BIGL(ep))) / (1-cr1);
end

function [ret] = X(x)
% at I1=L, envelop and cap intersect
% we have (X-L)/R = F(L), so X = L + R(L)*F(L)
DPCGlobals;
y = x + R(BIGL(x)) * F(x);
ret = max(0.0, y);
end

function [ret] = EVP(xl)
% one-to-one evp = W * (1 - exp(-D*X))
% if X=0, then evp=0
DPCGlobals;
ret = aw * (1 - exp(-ad*xl));
end

function [ret] = SJ2C(sj1,xl,capl)
% cap: 1/R * sqrt((X-L)^2 - (I1-L)^2)
DPCGlobals;
y2 = (xl-capl)^2 - (sj1-capl)^2;
rr = R(capl);
ret = sqrt(abs(y2)) / rr;
end

function [ret] = BMOD(sj1)
% elastic bulk
DPCGlobals;
ret = aki/(1-ak1) * (1 - ak1*exp(-ak2*sj1));
end

function [ret] = SMOD(sj2)
% elastic bulk
DPCGlobals;
ret = agi/(1-ag1) * (1 - ag1*exp(-ag2*sj2));
end





