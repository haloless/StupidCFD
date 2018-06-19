%
% must define bridge
%

%
% np = 41;
np = 61;
% np = 81;
[rp,xp] = AxisymEvolverGuessInit(bridge,np);

if 0
	% draw initial shape
	figure(hfig);
	hold on;
	plot(xp,rp,'-b');
	hold off;
	drawnow;
end
% return

% optimize
[rp,xp,pres] = AxisymEvolver(bridge,np,rp,xp);
%
pres

if 1
	% draw optimized shape
	figure(hfig);
	hold on;
	plot(xp,rp,'.-r');
	hold off;
	drawnow;
end 

% force
[F1,~,~,F2,~,~] = AxisymEvolverEvalForce(bridge, np,rp,xp,pres);



