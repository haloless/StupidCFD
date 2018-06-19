function [Dalg,sigma,epsE] = vmconst(epsEtr,E,v,fc)
% Constituive model for Prandtl-Reuss (von Mises perfect plasticity).
% So it has hardening H = 0
% use associate flow rule g = f
% 
% The system has unknowns by [eps,dgamma]
% The equations to be solved are
% 1~6: eps - eps_trial + dgamma*dg/dsigma = 0
% 7  : F(sigma) = 0
% 
% D is the consistent tangent modulus
% epsE is the elastic strain
% epsEtr: given trial elastic strain

tol = 1.0e-15;
maxit = 25;
bm1 = [1;1;1;0;0;0];

Ce = 1/E .* [-ones(3)*v+(1+v)*eye(3), zeros(3); zeros(3), 2*(1+v)*eye(3)];
De = inv(Ce);

% trial cauchy
sigma = De * epsEtr;
% deviatoric part
s = sigma - sum(sigma(1:3))/3*bm1;
% J2 invariant
j2 = (s'*s + s(4:6)'*s(4:6))/2;

% yield condition
f = sqrt(3*j2)/fc - 1;
% assume trial state
epsE = epsEtr;
Dalg = De;

if (f > tol)
	% yield, elasto-plastic
	% solve nonlinear material
	
	% residual vector, [eps; fyield]
	b = zeros(7,1);
	b(7) = f;
	
	% scalar plastic multiplier
	dgam = 0;
	
	dj2 = s;
	dj2(4:6) = 2*dj2(4:6);
	ddj2 = [eye(3)-ones(3)/3, zeros(3); zeros(3), 2*eye(3)];
	
	df = sqrt(3)/(2*fc*sqrt(j2)) * dj2;
	ddf = sqrt(3)/2/fc * (-dj2*dj2'/(2*j2^(3/2)) + ddj2/sqrt(j2));
	
	% material newton iteration
	itnum = 0;
	while (itnum<maxit) && ((norm(b(1:6))>tol) || (abs(b(7))>tol)
		itnum = itnum + 1;
		
		A = [eye(6)+dgam*ddf*De, df; df'*De, 0];
		dx = -A \ b;
		
		epsE = epsE + dx(1:6);
		dgam = dgam + dx(7);
		
		sigma = De * epsE;
		s = sigma - sum(sigma(1:3))/3*bm1;
		j2 = (s'*s + s(4:6)'*s(4:6))/2;
		
		dj2 = s;
		dj2(4:6) = 2*dj2(4:6);
		ddj2 = [eye(3)-ones(3)/3, zeros(3); zeros(3), 2*eye(3)];
		
		f = sqrt(3*j2)/fc - 1;
		df = sqrt(3)/(2*fc*sqrt(j2)) * dj2;
		ddf = sqrt(3)/2/fc * (-dj2*dj2'/(2*j2^(3/2)) + ddj2/sqrt(j2));
		
		% update residual
		b = [epsE-epsEtr+dgam*df; f];
	end
	
	% calc consistent tangent modulus
	B = inv([Ce, eye(6); dgam*ddf, -eye(6)]);
	Dalg = B(1:6,1:6) - (B(1:6,7:12)*(df*df')*B(1:6,1:6)) / ([df',zeros(1,6)]*B*[zeros(6,1); df]);
end





return
end

