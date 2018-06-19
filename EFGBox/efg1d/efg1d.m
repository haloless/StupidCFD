
clear;

%%
L = 1.0;

x = linspace(0,L,41).';
dx = x(2) - x(1);
xc = 0.5 * (x(1:end-1)+x(2:end));



nnodes = length(x);
ncells = nnodes - 1;

% 
E = 1.0;
area = 1.0;

%
dmax = 2.0;
dm = dmax * dx * ones(nnodes,1);

% gauss points, weights, jacobian for each cell
gg = [];
if 0
	jac = dx / 2;
	weight = 2;
	gg = xc;
end
if 0
	jac = dx / 2;
	weight = 1;
	gg = [xc-0.5*dx*sqrt(1/3); xc+0.5*dx*sqrt(1/3)];
end
if 1
    % use equidistance points for quadrature
    jac = dx / 2;
    nn = 2;
    weight = 2 / nn;
    for ii = 1:nn
        xx = xc - 0.5*dx + (ii-0.5)*(dx/nn);
        gg = [gg;xx];
    end
end

% add endpoint
gg = [0.0; gg];
gg = [gg;L];


if 1
	figure;
	plot(x,zeros(size(x)),'+-', gg,zeros(size(gg)),'x');
	pause;
end


%%
K = zeros(nnodes);
f = zeros(nnodes,1);
% G = zeros(nnodes,1);
G = zeros(nnodes,2);

for j = 1:length(gg)
	xg = gg(j);
	dif = xg - x;
	
	w = zeros(nnodes,1);
	dw = zeros(nnodes,1);
	
	for i = 1:nnodes
		drdx = sign(dif(i)) / dm(i);
		r = abs(dif(i)) / dm(i);
		if r <= 0.5
			w(i) = 2/3 - 4*r^2 + 4*r^3;
			dw(i) = (-8*r + 12*r^2) * drdx;
		elseif r <= 1.0
			w(i) = 4/3 - 4*r + 4*r^2 - 4/3*r^3;
			dw(i) = (-4 + 8*r - 4*r^2) * drdx;
		end
	end
	
	won = ones(1,nnodes);
	p = [ won; x.' ];
	B = p .* [w,w].';
	pp = zeros(2);
	A = zeros(2);
	dA = zeros(2);
	for i = 1:nnodes
		pp = p(:,i) * p(:,i)';
		A = A + w(i).*pp;
		dA = dA + dw(i).*pp;
	end
	Ainv = inv(A);
	
	pg = [1,xg];
	phi = pg * Ainv * B;
	db = p .* [dw,dw].';
	da = -Ainv * dA * Ainv;
	dphi = [0,1]*Ainv*B + pg*(da*B+Ainv*db);
	
	if j == 1 || j == length(gg)
		% endpoint, assemble G
		if j==1; jj=1; else; jj=2; end
		% G(1:3,jj) = -phi(1:3).';
		G(:,jj) = -phi.';
	else
		% integ point, assemble K,f
		K = K + (weight*E*area*jac)*(dphi'*dphi);
		% fbody = area * xg;
		fbody = area * 0;
		f = f + (weight*fbody*jac)*phi';
	end
end

% q = [0];
q = [0.0;-0.5];

mat = [ K, G; G', zeros(2) ];
rhs = [ f; q ];
sol = mat \ rhs;

u = sol(1:nnodes);

% uana = 1/E * (1/2*x - 1/6*x.^3);
uana = 0.5 * x;

if 1
	% figure;
	plot(x,uana,'-', x,u,'x');
	legend('ana','efg');
end

if 0
	% neig = 20;
	neig = nnodes - 2;
	% Ktan = mat;
	Ktan = K; Ktan(1,1) = 1.0e10;
	[v,d] = eigs(Ktan,neig,'sa');
	figure;
	for i = 1:neig
		vv = v(:,i);
		% vv = vv(1:end-1);
		dd = d(i,i);
		plot(x,vv);
		title(['eig(',int2str(i),')=',num2str(dd)]);
		pause;
	end
end





