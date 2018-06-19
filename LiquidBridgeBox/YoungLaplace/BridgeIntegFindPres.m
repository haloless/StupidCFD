function [res,pres] = BridgeIntegFindPres(bridge, alpha1,pres)
	
	R1 = bridge.R1;
	R2 = bridge.R2;
	H = bridge.H;
	theta1 = bridge.theta1;
	theta2 = bridge.theta2;
	X1 = bridge.X1;
	X2 = bridge.X2;
	
	pres = nan;
	plo = nan;
	phi = nan;
	elo = 0;
	ehi = 0;
	
	%
	for pp = 1.3:0.001:1.5
		rr = BridgeInteg(bridge, alpha1, pp);
		if rr.touch
			ee = rr.ca2 - theta2;
			[pp,ee]
			if isnan(plo)
				plo = pp;
				elo = ee;
			elseif ee*elo <= 0
				phi = pp;
				ehi = ee;
				break;
			% elseif ee*elo > 0 && abs(ee)<abs(elo)
				% plo = pp;
				% elo = ee;
			end
		else
			if ~isnan(plo) && ~isnan(phi)
				break;
			end
		end
	end
	
	plo
	phi
	
	if isnan(plo) || isnan(phi)
		return
	end
	
	res = [];
	errmin = 1.0e10;
	for pp = plo:1e-2:2.0
		rr = BridgeInteg(bridge, alpha1,pp);
		if (rr.touch)
			caerr = abs(rr.ca2-theta2);
			if caerr < errmin
				errmin = caerr;
				res = rr;
				pres = pp;
			end
		end
	end
	


return
end



