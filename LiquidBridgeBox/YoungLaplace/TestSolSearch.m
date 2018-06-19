
function [] = TestSolSearch(R1,R2,theta1,theta2, V, Hbegin, varargin)
    
    sigma = 1.0;
    
    theta1 = deg2rad(theta1);
    theta2 = deg2rad(theta2);
    
    H = Hbegin;
    
    % initial bridge
    bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
    
    % draw spheres
    hfig = figure;
    hold on;
    PlotBridgeGeom(bridge);
    hold off;
    
    % optimize
    np = 51;
    % np = 81;
    % np = 101;
    [rp,xp] = AxisymEvolverGuessInit(bridge,np);
    [rp,xp,pres] = AxisymEvolver(bridge,np,rp,xp);
    
    if 1
        % draw optimized shape
        figure(hfig);
        hold on;
        plot(xp,rp,'.-r');
        hold off;
        drawnow;
    end 
    
    % initial state
    alpha1 = asin(rp(1)/R1);
    if R2 > 0
        alpha2 = asin(rp(end)/R2);
    else
        alpha2 = rp(end);
    end
    pres = pres;
    
    if 0
        % initial integral solution 
        sol = BridgeInteg(bridge, alpha1,pres);
    else
        sol = BridgeIntegFindAlfa(bridge, pres, alpha1,alpha2,H);
        alpha1 = sol.alpha1;
        alpha2 = sol.alpha2;
        H = sol.H;
    end
    
    if 1
        figure(hfig);
        hold on;
        plot(sol.xs,sol.ys,'-b');
        hold off;
        drawnow;
    end
    if 1
        pause;
    end
    
    %
    data = [];
    
    Hbeg = H;
    Hmax = H;
    solmax = [];
    brgmax = [];
    
    %
    pincr = 1e-2;
    for i = 1:2:length(varargin)
        key = varargin{i};
        val = varargin{i+1};
        switch key
        case 'pincr'
            pincr = val;
        otherwise
            error('Unknown %s=%s', key,val);
        end
    end
    
    
    for iter = 1:100000
        
        
        pnew = pres + pincr;
        a1old = alpha1;
        a2old = alpha2;
        Hold = H;
        
        %
        [soltr] = BridgeIntegFindAlfa(bridge, pnew, a1old,a2old,Hold);
        
        %
        sol = soltr;
        
        pres = sol.pres;
        alpha1 = sol.alpha1;
        alpha2 = sol.alpha2;
        H = sol.H;
        
        %
        bridge = MakeBridge(R1,R2,H,theta1,theta2,V,sigma);
        
        % save data
        data(end+1,:) = [H,alpha1,pres];
        
        if H > Hmax
            Hmax = H;
            solmax = sol;
            brgmax = bridge;
        end
        
        if mod(iter,10) == 0 
            disp(sprintf('iter=%d, pres=%0.4f, H=%0.4f', iter,pres,H));
            
            figure(hfig);
            plot(sol.xs,sol.ys,'-b');
            hold on;
            PlotBridgeGeom(bridge);
            hold off;
            title(sprintf('P=%0.6f,H=%0.6f', pres, H));
            drawnow;
        end
        
        if H<Hmax - 0.05*(Hmax-Hbeg)
            break;
        end
    end
    
    %
    [Hmax,ind] = max(data(:,1));
    % data(ind-1:ind+1,1)
    disp(['Hmax=',num2str(Hmax)]);
    a1max = data(ind,2);
    disp(['a1max=',num2str(a1max)]);
    b1max = sin(a1max);
    disp(['b1max=',num2str(b1max)]);
    if R2 <= 0
        disp(['b2max=',num2str(solmax.alpha2)]);
    end
    
    solmax
    (solmax.xneck-R1)/solmax.H
    
    if 1
        % alpha1 vs H
        figure;
        subplot(1,2,1);
        plot(data(:,1),data(:,2),'.-', data(ind,1),data(ind,2),'rx');
        % vs P
        % figure;
        subplot(1,2,2);
        plot(data(:,3),data(:,1),'.-', data(:,3),data(:,2),'.-', data(ind,3),data(ind,1),'x', data(ind,3),data(ind,2),'x');
        legend('H','a1');
    end
    if 1
        figure(hfig);
        plot(solmax.xs,solmax.ys,'-b', solmax.xe,solmax.ye,'xr','MarkerSize',10);
        hold on;
        PlotBridgeGeom(brgmax);
        hold off;
    end
return
end

