function [ts,zs, res,exitflag] = YLSolvePres2(bridge, alpha1, pguess)
%YLSolvePres2 
    
    R1 = bridge.R1;
    R2 = bridge.R2;
    H = bridge.H;
    theta1 = bridge.theta1;
    theta2 = bridge.theta2;
    
    %
    function [eout] = bridgefunc(pin)
        
        eout = zeros(size(pin));
        
        for i = 1:numel(pin)
            pp = pin(i);
            [ts,zs, res] = YLIntegrateBridge(bridge, alpha1, pp);
            x2 = ts(end);
            y2 = zs(end,1);
            dy2 = zs(end,2);
            
            ee = sqrt((x2-res.x2)^2 + (y2-res.y2)^2 + (dy2-res.dy2)^2);
            
            eout(i) = ee;
        end
        
    return
    end
    
    %
    options = optimset('Display','off');
    [pres,residual,exitflag] = fsolve(@bridgefunc, pguess, options);
    
    [ts,zs, res] = YLIntegrateBridge(bridge, alpha1,pres);
    
    res.pres = pres;



return
end







