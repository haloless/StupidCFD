%testTubeHillElastoPlastic: thick cylinder problem
% Hill's elastoplastic solution

Y = Fyield * 2 / sqrt(3);
% pressure begin to yield
Pinit = Y/2 * (1-Ra^2/Rb^2);
% pressure fully yield
Pfinal = Y * log(Rb/Ra);
if P0 > Pfinal
    % fully yield
    Rc = Rb;
    disp('Plastic');
elseif P0 > Pinit
    % partial yield radius
    Rcfun = @(x) log(x/Ra) + 0.5*(1-x.^2/Rb^2) - P0/Y;
    Rc = fsolve(Rcfun, (Ra+Rb)/2);
    disp(['Elastoplastic, Rc=',num2str(Rc)]);
else
    Rc = Ra;
    disp('Elastic');
end

rcheck = [Ra:0.005:Rc,Rc:0.005:Rb];
sigradial = zeros(size(rcheck));
sighoop = zeros(size(rcheck));
sigradial(rcheck<Rc) = Y.*(-0.5-log(Rc./rcheck(rcheck<Rc))+0.5*Rc^2/Rb^2);
sighoop(rcheck<Rc) = Y.*(0.5-log(Rc./rcheck(rcheck<Rc))+0.5*Rc^2/Rb^2);
sigradial(rcheck>=Rc) = -0.5*Y*Rc^2/Rb^2 .* (Rb^2 ./ rcheck(rcheck>=Rc).^2 - 1);
sighoop(rcheck>=Rc) = 0.5*Y*Rc^2/Rb^2 .* (Rb^2 ./ rcheck(rcheck>=Rc).^2 + 1);

if 1 % plot elastic/plastic nodes
    figure;
    epmask = (nodeFlag==1);
    plot(nodeCurr(1,epmask),nodeCurr(2,epmask),'x', nodeCurr(1,~epmask),nodeCurr(2,~epmask),'o')
    axis equal;
    legend('plastic','elastic');
end

if 1 % plot hoop stress
    figure;
    if prob_type == 1
        % (syy on X-axis)
        iplot = find(nodeX>0 & abs(nodeY)<1.0e-6);
        plot(rcheck,sighoop,'-', nodeX(iplot),nodeSigma(2,iplot),'x');
    elseif prob_type == 2
        % (szz in RZ coord)
        iplot = find(nodeX>0 & abs(nodeY-H/2)<h0);
        plot(rcheck,sighoop,'-', nodeX(iplot),nodeSigma(4,iplot),'x');
    end
    title('hoop stress');
end

if 1 % plot radial stress
    figure;
    if prob_type == 1
        % (sxx on X-axis)
        iplot = find(nodeX>0 & abs(nodeY)<1e-6);
        plot(rcheck,sigradial,'-', nodeX(iplot),nodeSigma(1,iplot),'x');
    elseif prob_type == 2
        % (sxx in RZ coord)
        iplot = find(nodeX>0 & abs(nodeY-H/2)<h0);
        plot(rcheck,sigradial,'-', nodeX(iplot),nodeSigma(1,iplot),'x');
    end
    title('radial stress')
end


