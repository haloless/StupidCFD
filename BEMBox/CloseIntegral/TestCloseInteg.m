
clear all;

L = 0.5;

% points = [...
% -L, -L;
 % L, -L;
 % L,  L;
% -L,  L];

x0 = [ 0; 0; ];
x1 = [ -L; -L; ];
x2 = [  L; -L; ];
x3 = [  L;  L; ];
x4 = [ -L;  L; ];

qfunc = @(x,y) 1.0 ./ sqrt(x.^2+y.^2);

% matlab quadrature
qref = integral2(qfunc, -L,L, -L,L)

ng = 7;
[wg,xg,yg] = TriangleGaussRule(ng);

% use direct gauss quadrature
qint = 0;
qint = qint + TriangleGaussInteg(qfunc,x0,x1,x2, ng,wg,xg,yg);
qint = qint + TriangleGaussInteg(qfunc,x0,x2,x3, ng,wg,xg,yg);
qint = qint + TriangleGaussInteg(qfunc,x0,x3,x4, ng,wg,xg,yg);
qint = qint + TriangleGaussInteg(qfunc,x0,x4,x1, ng,wg,xg,yg);

%
ng = 7;
[wg,xg] = LineGaussRule(ng);

qrad = 0;
qrad = qrad + TriangleRadSingularInteg(qfunc,x0,x1,x2, ng,wg,xg, ng,wg,xg);
qrad = qrad + TriangleRadSingularInteg(qfunc,x0,x2,x3, ng,wg,xg, ng,wg,xg);
qrad = qrad + TriangleRadSingularInteg(qfunc,x0,x3,x4, ng,wg,xg, ng,wg,xg);
qrad = qrad + TriangleRadSingularInteg(qfunc,x0,x4,x1, ng,wg,xg, ng,wg,xg);

