
if 1
    
    % figure(hfig);
    % hold on;
    
    % plot cap
    a = par.a0 + alpha*par.Hcap;
    plotCap(par, a, 'r-');
    
    % plot current state
    [p,s] = voigt3dPressShear(sigma);
    sj2 = sqrt(voigt3dJ2(s));
    plot(p,sj2,'x');
    
    titlestr = ['step=',int2str(istep),'/',int2str(nstep)];
    switch mtype
    case 0
        titlestr = [titlestr,';mtype=',int2str(mtype),'(elastic)'];
    case 1
        titlestr = [titlestr,';mtype=',int2str(mtype),'(cone)'];
    case 2
        titlestr = [titlestr,';mtype=',int2str(mtype),'(apex)'];
    case 3
        titlestr = [titlestr,';mtype=',int2str(mtype),'(cap)'];
    case 4
        titlestr = [titlestr,';mtype=',int2str(mtype),'(corner)'];
    end
    title(titlestr);
    
    % hold off;

end

