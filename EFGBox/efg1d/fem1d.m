
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
ng = 2;
gg = [];
if 0
	jac = dx / 2;
	weight = 2;
	gg = xc;
end
if ng == 2
	jac = dx / 2;
	weight = 1;
	gg = [xc-0.5*dx*sqrt(1/3), xc+0.5*dx*sqrt(1/3)];
end
if 0
    jac = dx / 2;
    nn = 50;
    weight = 2 / nn;
    for ii = 1:nn
        xx = xc - 0.5*dx + (ii-0.5)*(dx/nn);
        gg = [gg;xx];
    end
end

% add endpoint
nlag = 2;
% gg = [0.0; gg];
if nlag > 1
    % gg = [gg;L];
end


if 1
	figure;
	plot(x,zeros(size(x)),'+-', gg(:),zeros(size(gg(:))),'x');
	pause;
end


%%
K = zeros(nnodes);
f = zeros(nnodes,1);
G = zeros(nnodes,2);

for i = 1:ncells
    for j = 1:ng
        xg = gg(i,j);
        
        ia = i;
        ib = i+1;
        xa = x(ia);
        xb = x(ib);
        
        phi = zeros(nnodes,1);
        phi(ia) = (xb-xg)/dx;
        phi(ib) = (xg-xa)/dx;
        
        dphi = zeros(nnodes,1);
        dphi(ia) = -1/dx;
        dphi(ib) = 1/dx;
        
        % integ point, assemble K,f
        K = K + (weight*E*area*jac)*(dphi*dphi');
        % fbody = area * xg;
        fbody = area * 0;
        f = f + (weight*fbody*jac)*phi;
    end
end

G(1,1) = -1;
G(nnodes,2) = -1;


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

if 1
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





