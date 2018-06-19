% Plot radial displacement
% NOTE this is the elastic solution


uv = reshape(uv, 2,numNodes);

nodeCurr = nodeCoord + uv;

% current config
figure;
plot(nodePos(:,1),nodePos(:,2),'o', nodeCurr(1,:),nodeCurr(2,:),'x');
legend('before','after');
axis equal;

% radial displacement
figure;
if prob_type == 1
    rr2 = sqrt(nodeX.^2+nodeY.^2)';
    ur2 = (uv(1,:).*nodeX' + uv(2,:).*nodeY')./rr2;
elseif prob_type == 2
    rr2 = nodeX.';
    ur2 = uv(1,:);
end
for irad = 1:nrad
    rr(irad) = Ra + (irad-0.5)*h0;
    ur(irad) = mean(ur2(abs(rr2-rr(irad))<1.0e-6));
end

% analytical plane-strain solution
c1 = (1+nu0)/E0*(1-2*nu0)*P0/(Rb^2/Ra^2-1);
c2 = (1+nu0)/E0*Rb^2*P0/(Rb^2/Ra^2-1);
rs = linspace(Ra,Rb,100);
us = c1 .* rs + c2 ./ rs;
% plot(rr,ur,'x',rs,us,'-');
plot(rr2,ur2,'x',rs,us,'-');
legend('sim','ana');
title('Radial disp.');

