
%% solve infinitesimal problem
if 1
    fext = zeros(numDofs,1);
    fext(tracBCDofs) = tracBCVals;
    
    frhs = fext - fint;
    
    sol = fdpmSolve(Ktan0,frhs, dispBCDofs,dispBCVals);
    
    uv = reshape(sol, 2,numNodes);
    
    % update particle coordinates
    nodeCoord = nodeCoord + uv * 10;
    % xpos = reshape(nodeCoord(1,:), nx,ny);
    % ypos = reshape(nodeCoord(2,:), nx,ny);
    
    % plot particles after displacement
    figure;
    plot(nodePos(:,1),nodePos(:,2),'o', nodeCoord(1,:),nodeCoord(2,:),'x');
    legend('before','after');
    axis equal;
    
    
    % plot radial displacement
    figure;
    if prob_type == 1
        % plane-strain config
        % average radial disp for each layer
        rr2 = sqrt(nodeX.^2+nodeY.^2)';
        ur2 = (uv(1,:).*nodeX' + uv(2,:).*nodeY')./rr2;
        for irad = 1:nrad
            rr(irad) = Ra + (irad-0.5)*h0;
            ur(irad) = mean(ur2(abs(rr2-rr(irad))<1.0e-6));
        end
    elseif prob_type == 2
        % axisymmetric config
        rr2 = nodeX.';
        ur2 = uv(1,:);
        for irad = 1:nrad
            rr(irad) = Ra + (irad-0.5)*h0;
            ur(irad) = mean(ur2(abs(rr2-rr(irad))<1.0e-6));
        end
    end
    
    % analytical plane-strain solution
    c1 = (1+nu0)/E0*(1-2*nu0)*P0/(Rb^2/Ra^2-1);
    c2 = (1+nu0)/E0*Rb^2*P0/(Rb^2/Ra^2-1);
    rs = linspace(Ra,Rb,100);
    us = c1 .* rs + c2 ./ rs;
    plot(rr,ur,'x',rs,us,'-');
    
end

% return
