function [L] = parDerGen(X,eV,eP,yP,yd)
%parDerGen: Partial Derivative of a General 2nd-order tensor function w.r.t its argument
% Return tensor L(ij)(kl)
% Inputs of parDerGen
% X is the 3x3 tensor argument
% eV is its eigenvector 3x3
% eP is its eigenvalue 3x1
% yP is the function applied to eigenvalues
% yd is the function derivative to eigenvalues
% Outputs of parDerGen
% L is 9x9 tensor
% The order is [ xx,yy,zz, xy,yx, yz,zy, zx,xz]
%



tol = 1.0e-9;

s = zeros(5,1);
Is = [eye(3), zeros(3); zeros(3), eye(3)/2];
bm1 = [1,1,1,0,0,0]';

if abs(eP(1))<tol && abs(eP(2))<tol && abs(eP(3))<tol
	% eigenvalues close to zero
	L = Is;
elseif abs(eP(1)-eP(2))<tol && abs(eP(1)-eP(3))<tol
	% 3 equal eigenvalues
	L = yd(1) .* Is;
elseif abs(eP(1)-eP(2))<tol || abs(eP(2)-eP(3))<tol || abs(eP(1)-eP(3))<tol
	% 2 equal eigenvalues
	x = [X(1),X(5),X(9),X(2),X(6),X(7)]';
	
	s(1) = (yP(1)-yP(3)) / (eP(1)-eP(3))^2 - yd(3)/(eP(1)-eP(3));
	s(2) = 2*eP(3)*(yP(1)-yP(3))/(eP(1)-eP(3))^2 - (eP(1)+eP(3))/(eP(1)-eP(3))*yd(3);
	s(3) = 2*(yP(1)-yP(3))/(eP(1)-eP(3))^3 - (yd(1)+yd(3))/(eP(1)-eP(3))^2;
	s(4) = eP(3) * s(3);
	s(5) = eP(3)^2 * s(3);
	
	dX2dX = [...
	2*X(1), 0, 0, X(2), 0, X(3);
	0, 2*X(5), 0, X(2), X(6), 0;
	0, 0, 2*X(9), 0, X(6), X(3);
	X(2), X(2), 0, (X(1)+X(5))/2, X(3)/2, X(6)/2;
	0, X(6), X(6), X(3)/2, (X(5)+X(9))/2, X(2)/2;
	X(3), 0, X(3), X(6)/2, X(2)/2, (X(1)+X(9))/2;
	];
	
	L = s(1)*dX2dX - s(2)*Is - s(3)*(x*x') + s(4)*(x*bm1'+bm1*x') - s(5)*(bm1*bm1');
else
	% all different eigenvalues
	D = [(eP(1)-eP(2))*(eP(1)-eP(3)); (eP(2)-eP(1))*(eP(2)-eP(3)); (eP(3)-eP(1))*(eP(3)-eP(2))];
	alfa = 0;
	bta = 0;
	g = zeros(3,1);
	eD = zeros(6,3);
	for i = 1:3
		alfa = alfa + yP(i)*eP(i)/D(i);
		bta = bta + yP(i)/D(i)*det(X);
		for j = 1:3
			g(i) = g(i) + yP(j)*eP(j)/D(j)*(det(X)/eP(j)-eP(i)^2)/eP(i)^2;
		end
		esq = eV(:,i)*eV(:,i)';
		eD(:,i) = [esq(1,1),esq(2,2),esq(3,3),esq(1,2),esq(2,3),esq(3,1)]';
	end
	
	y = inv(X);
	
	Ib = [...
	y(1)^2, y(2)^2, y(7)^2, y(1)*y(2), y(2)*y(7), y(1)*y(7);
	y(2)^2, y(5)^2, y(6)^2, y(5)*y(2), y(5)*y(6), y(2)*y(6);
	y(7)^2, y(6)^2, y(9)^2, y(6)*y(7), y(9)*y(6), y(9)*y(7);
	y(1)*y(2), y(5)*y(2), y(6)*y(7), (y(1)*y(5)+y(2)^2)/2,    (y(2)*y(6)+y(5)*y(7))/2, (y(1)*y(6)+y(2)*y(7))/2;
	y(2)*y(7), y(5)*y(6), y(9)*y(6), (y(2)*y(6)+y(5)*y(7))/2, (y(9)*y(5)+y(6)^2)/2,    (y(9)*y(2)+y(6)*y(7))/2;
	y(1)*y(7), y(2)*y(6), y(9)*y(7), (y(1)*y(6)+y(2)*y(7))/2, (y(9)*y(2)+y(6)*y(7))/2, (y(9)*y(1)+y(7)^2)/2;
	];
	
	L = alfa*Is - bta*Ib;
	L = L + (yd(1)+g(1))*eD(:,1)*eD(:,1)' + (yd(2)+g(2))*eD(:,2)*eD(:,2)' + (yd(3)+g(3))*eD(:,3)*eD(:,3)';
	% L = L + eD*diag(yd+g)*eD';
end

L = [L(1:3,1:3), L(1:3,[4,4,5,5,6,6]); L([4,4,5,5,6,6],1:3), L([4,4,5,5,6,6],[4,4,5,5,6,6])];

return
end


