
if (0)
	% check
	for i = 1:numNodes
		tx = eye(2,2);
		u = nodePos * tx';
		Fx = u' * [dNX(i,:)',dNY(i,:)'];
		if (norm(Fx-tx) > 1e-8)
			disp(['check1: node ', int2str(i), ' dXdX error']);
			disp([Fx tx]);
		end
		
		% nodex = nodeX * 1 + nodeY * 2;
		% nodey = nodeX * 3 + nodeY * 4;
		tx = [1 2;3 4];
		u = nodePos * tx';
		Fx = u' * [dNX(i,:)',dNY(i,:)'];
		if (norm(Fx-tx) > 1e-8)
			disp(['check2: node ', int2str(i), ' dxdX error']);
			disp([Fx tx]);
		end
	end
end

if (0)
	% this is for order(2) only
	for i = 1:numNodes
		u = nodePos.^2;
		Fx = u' * [dNX(i,:)',dNY(i,:)'];
		tx = [2*nodePos(i,1),0; 0,2*nodePos(i,2)];
		if (norm(Fx-tx) > 1e-8)
			disp(['check3: node ', int2str(i), ' dxdX error']);
			disp([Fx tx]);
		end
		
		u = [1/2*nodePos(:,1).^2, nodePos(:,1).*nodePos(:,2), 1/2*nodePos(:,2).^2];
		Fx = u' * [dNXX(i,:)', dNXY(i,:)', dNYY(i,:)'];
		tx = [1 0 0; 0 1 0; 0 0 1];
		if (norm(Fx-tx) > 1e-8)
			disp(['check4: node ', int2str(i), ' dxdX error']);
			disp([Fx tx]);
		end
	end
end


