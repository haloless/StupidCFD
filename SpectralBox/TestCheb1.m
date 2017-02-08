
clear all;

fun = @(x) 1.0 ./ (1+16*x.^2);

% analytical
xx = -1.01:.005:1.01;
uu = fun(xx);

% N = 8;
% N = 16;
% N = 24;
% N = 32;
% N = 64;
N = 128;
% N = 192;
% N = 256;

% for i = 1:2
% if i==1, s = 'equispaced points'; x = -1 + 2*(0:N)/N; end
% if i==2
if 1
	s = 'Chebyshev points'; 
	x = cos(pi*(0:N)/N); 
end

% collocation value
u = fun(x);

if 0
	p = polyfit(x,u,N); % interpolation
	pp = polyval(p,xx); % evaluation of interpolant
end

if 0
	A = chebpoly1(N,x');
	c = A \ u';
	
	aa = chebpoly1(N,xx');
	pp = aa * c;
	pp = pp';
end

if 1
	A = chebpoly1(N,x');
	
	% c_hat = 2 for k=0,N; =1 for 0<k<N
	% this is the corresponding coefficients of discrete orthogonality
	% for Chebyshev 1st-polynomial on Chebyshev 2nd-node
	chat = [ 2; ones(N-1,1); 2 ];
	
	uhat = u' ./ chat;
	
	ahat = A' * uhat;
	
	c = 2/N * ahat ./ chat;
	
	aa = chebpoly1(N,xx');
	pp = aa * c;
	pp = pp';
end

figure;

% plot(x,u,'.','markersize',13)
% line(xx,pp,'linewidth',.8)
plot(x,u,'.', xx,pp,'-', xx,uu,'--', 'markersize',16); 
legend('point','chebpoly','fun');
axis([-1.1 1.1 -1 1.5])
title(s)
error = norm(uu-pp,inf);
text(-0.5,-0.5,['max error = ' num2str(error)])

