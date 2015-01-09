
% clc
clear all

% ader_order = 3;
ader_order = 4;
% ader_order = 5;

ADERWENOGlobals1D;
ADERWENOInit1D(ader_order);


% M = 2; % 3rd-order
% M = 3; % 4th-order
% M = 4; % 5th-order
% N = M + 1;
M = MDegree;
N = NPoint;

eta = GausEta;

% % number of stencils
% if (mod(N,2) == 1)
    % Ns = 3; 
% else
    % Ns = 4;
% end

% % Gaussian quadrature points and weights on [0,1] range
% [eta,wgt] = GaussQuadCoefs1D(N, 0, 1);

% % stencil offset
% if (mod(N,2) == 1) 
    % soff = [-M, -M/2, 0];
% else
    % soff = [-M, -N/2, -N/2+1, 0];
% end
% % coefficient matrix of reconstruction
% smat = zeros(N,N,Ns);

% lagpoly = LagInterpPoly(eta);

% for s = 1:Ns
    % cmat = zeros(N,N);
    % for e = 0:M
        % xe = eta + e+soff(s);
        % ce = zeros(N,N);
        % for i = 1:N
            % % ce(i,:) = LagInterpCoef(xe(i), eta');
            % for j = 1:N
                % ce(i,j) = polyval(lagpoly(:,j), xe(i));
            % end
        % end
        % cmat(e+1,:) = wgt' * ce;
    % end
    % smat(:,:,s) = cmat;
% end

% % oscillation indicator matrix
% Sigmat = zeros(N,N);
% derlagp = lagpoly;
% for a = 1:M
    % derlagp1 = zeros(N-a,N);
    % for p = 1:N
        % derlagp1(:,p) = polyder(derlagp(:,p));
    % end
    % derlagp = derlagp1;
    % for p = 1:N
    % for m = 1:N
        % vapm = sum(polyval(derlagp(:,p),eta) .* polyval(derlagp(:,m),eta) .* wgt);
        % Sigmat(p,m) = Sigmat(p,m) + vapm;
    % end
    % end
% end

% test
% f = @(x) (2*x+1);
% xcen = 0.5;
% hx = 0.2;

% f = @(x) (x.^2 + 2*x + 3);
% xcen = -1.0;
% hx = 0.01;

% f = @(x) (10*x.^2 + 2*x + 3);
% xcen = -1.0;
% hx = 0.01;

% f = @(x) (-10*x.^3 + x.^2 - 2*x + 3);
% xcen = 0.0;
% hx = 0.01;

% discontinuous function
xcen = 0.0;
hx = 0.01;
% xstep = xcen;
% xstep = xcen + hx*0.5;
xstep = xcen - hx*0.5;
f = @(x) ((x>xstep).*1 + (1-(x>xstep)).*0);

icen = M+1;
xs = linspace(xcen-M*hx,xcen+M*hx,M*2+1);
if (0)
    % initialize as cell point value
    fs = f(xs);
else
    % initialize as cell average
    for i = icen-M:icen+M
        fs(i) = 1/hx * integral(f, xs(i)-hx/2, xs(i)+hx/2);
    end
end

qx = (xs(icen)-hx*0.5) + eta*hx;

% qs = zeros(N,Ns);
% ps = zeros(N,Ns);
% for s = 1:Ns
    % sidx = icen + (0:M);
    % sidx = sidx + soff(s);
    % rhs = fs(sidx)';
    % % computed point values
    % qs(:,s) = smat(:,:,s) \ rhs;
    
    % % Lagrange interpolating polynomials
    % if (1)
        % % ps(:,s) = lagpoly * qs(:,s);
        % ps(:,s) = SumPoly(lagpoly, qs(:,s));
    % else
        % ps(:,s) = LagInterpPoly(qx,qs(:,s));
    % end
% end


% % stencil weights
% sigs = zeros(Ns,1);
% for s = 1:Ns
    % sigs(s) = qs(:,s)' * Sigmat * qs(:,s);
% end
% if (mod(N,2) == 1)
    % lambdas = [1; 1e5; 1];
% else
    % lambdas = [1; 1e5; 1e5; 1];
% end
% epsil = 1e-14; wr = 8;
% omegas = lambdas ./ ((sigs+epsil).^wr);
% ws = omegas ./ sum(omegas);

% % WENO polynomial
% qweno = qs * ws;
% pweno = SumPoly(ps,ws);

[qweno,pweno,qs,ps,sigmas,omegas,ws] = ADERWENOReconstruct1D(fs(icen-M:icen+M),hx);

if (1)
    figure;
    xs1 = linspace(xs(1),xs(end),32);
    qeta1 = (xs1-(xcen-hx*0.5)) ./ hx;
    if (NStencil == 3) 
        plot(xs1,f(xs1),'-', xs,fs,'.', ...
            qx,qs(:,1),'x', qx,qs(:,2),'o', qx,qs(:,3),'s', ...
            qx,qweno,'d', xs1,polyval(pweno,qeta1),'--');
        legend('f','fbar','left','center','right', ...
            'weno-node','weno-poly');
    elseif (NStencil == 4)
        plot(xs1,f(xs1),'-', xs,fs,'.', ...
        qx,qs(:,1),'x', qx,qs(:,2),'o', qx,qs(:,3),'s', qx,qs(:,4),'^', ...
        qx,qweno,'d', xs1,polyval(pweno,qeta1),'--');
        legend('f','fbar','left','center_{left}','center_{right}','right', ...
        'weno-node','weno-poly');
    end
    
    for s = 1:3
        qerr(s) = max(abs(qs(:,s) - f(qx)));
    end
end


