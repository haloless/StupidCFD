function waterwave
% WATERWAVE   2D Shallow Water Model
%
% Lax-Wendroff finite difference method.
% Reflective boundary conditions.
% Random water drops initiate gravity waves.
% Surface plot displays height colored by momentum.
% Plot title shows t = simulated time and tv = a measure of total variation.
% An exact solution to the conservation law would have constant tv.
% Lax-Wendroff produces nonphysical oscillations and increasing tv.
%
% See:
%    http://en.wikipedia.org/wiki/Shallow_water_equations
%    http://www.amath.washington.edu/~rjl/research/tsunamis
%    http://www.amath.washington.edu/~dgeorge/tsunamimodeling.html
%    http://www.amath.washington.edu/~claw/applications/shallow/www

% Parameters

Lx = 256;
Ly = 256;
% n = 64;                  % grid size
n = 128;                  % grid size
% n = 256;                  % grid size
g = 9.8;                 % gravitational constant
dt = 0.01;               % hardwired timestep
% dx = 1.0;
% dy = 1.0;
dx = Lx / n;
dy = Ly / n;
nplotstep = 8;           % plot interval
ndrops = 1;              % maximum number of drops
dropstep = 200;          % drop interval
% D = droplet(1.5,21);     % simulate a water drop
D = droplet(2.5,21);     % simulate a water drop

Xs = linspace(-dx/2,dx*n+dx/2,n+2);
Ys = linspace(-dy/2,dy*n+dy/2,n+2);
[Xs,Ys] = ndgrid(Xs,Ys);

% Initialize graphics

[surfplot,top,restart,quit,bottomplot] = initgraphics(n);

% Outer loop, restarts.

while get(quit,'value') == 0
   set(restart,'value',0)
   
   H = ones(n+2,n+2);   U = zeros(n+2,n+2);  V = zeros(n+2,n+2);
   % Hx = zeros(n+1,n+1); Ux = zeros(n+1,n+1); Vx = zeros(n+1,n+1);
   % Hy = zeros(n+1,n+1); Uy = zeros(n+1,n+1); Vy = zeros(n+1,n+1);
   Hx = zeros(n+2,n+2); Ux = zeros(n+2,n+2); Vx = zeros(n+2,n+2);
   Hy = zeros(n+2,n+2); Uy = zeros(n+2,n+2); Vy = zeros(n+2,n+2);
   Hnew = zeros(n+2,n+2);
   Unew = zeros(n+2,n+2);
   Vnew = zeros(n+2,n+2);
   Fcx = zeros(n+2,n+2);
   Fcy = zeros(n+2,n+2);
   
   B = zeros(n+2,n+2); % bed height
   if (0)
       B(:,:) = 0.01*Xs(:,:) + 0.02*Ys(:,:);
   end
   B(:,1) = B(:,2); B(:,n+2) = B(:,n+1);
   B(1,:) = B(2,:); B(n+2,:) = B(n+1,:);
   S = zeros(n+2,n+2); % source term
   
   ndrop = ceil(rand*ndrops);
   nstep = 0;

   % Inner loop, time steps.

    while get(restart,'value')==0 && get(quit,'value')==0
        nstep = nstep + 1;

        % Random water drops
        if mod(nstep,dropstep) == 0 && nstep <= ndrop*dropstep
           w = size(D,1);
           i = ceil(rand*(n-w))+(1:w);
           j = ceil(rand*(n-w))+(1:w);
           H(i,j) = H(i,j) + (1+4*rand)/5*D;
        end

        % Reflective boundary conditions
        H(:,1) = H(:,2);      U(:,1) = U(:,2);       V(:,1) = -V(:,2);
        H(:,n+2) = H(:,n+1);  U(:,n+2) = U(:,n+1);   V(:,n+2) = -V(:,n+1);
        H(1,:) = H(2,:);      U(1,:) = -U(2,:);      V(1,:) = V(2,:);
        H(n+2,:) = H(n+1,:);  U(n+2,:) = -U(n+1,:);  V(n+2,:) = V(n+1,:);
        
        I = 2:n+1;
        J = 2:n+1;
        Is = 1:n+1;
        Js = 1:n+1;
        
        if (0)
            % Lax-Friedrichs
            % a = 0.995; % good for B=const
            a = 0.8;
            
            Hnew(I,J) = a * H(I,J) ...
            + (1-a)/4 * (H(I+1,J) + H(I-1,J) + H(I,J+1) + H(I,J-1)) ...
            - dt/(2*dx) * (U(I+1,J) - U(I-1,J)) ...
            - dt/(2*dy) * (V(I,J+1) - V(I,J-1));
            
            Unew(I,J) = a * U(I,J) ...
            + (1-a)/4 * (U(I+1,J) + U(I-1,J) + U(I,J+1) + U(I,J-1)) ...
            - dt/(2*dx) * ((U(I+1,J).^2)./H(I+1,J) + g/2*H(I+1,J).^2) ...
            + dt/(2*dx) * ((U(I-1,J).^2)./H(I-1,J) + g/2*H(I-1,J).^2) ...
            - dt/(2*dy) * (U(I,J+1).*V(I,J+1)./H(I,J+1) - U(I,J-1).*V(I,J-1)./H(I,J-1)) ...
            + dt * (-g/(2*dx)*H(I,J).*(B(I+1,J)-B(I-1,J)));
            
            Vnew(I,J) = a * V(I,J) ...
            + (1-a)/4 * (V(I+1,J) + V(I-1,J) + V(I,J+1) + V(I,J-1)) ...
            - dt/(2*dx) * (U(I+1,J).*V(I+1,J)./H(I+1,J) - U(I-1,J).*V(I-1,J)./H(I-1,J)) ...
            - dt/(2*dy) * ((V(I,J+1).^2)./H(I,J+1) + g/2*H(I,J+1).^2) ...
            + dt/(2*dy) * ((V(I,J-1).^2)./H(I,J-1) + g/2*H(I,J-1).^2) ...
            + dt * (-g/(2*dy)*H(I,J).*(B(I,J+1)-B(I,J-1)));
        end % LF
        
        if (1) % Lax-Wendroff
            % First half step
            % x direction
            i = 1:n+1;
            j = 1:n;
            % height
            Hx = (H(i+1,j+1)+H(i,j+1))/2 - dt/(2*dx)*(U(i+1,j+1)-U(i,j+1));
            % x momentum
            Fx = (U.^2)./H + g/2*H.^2;
            Ux = 1/2*(U(i+1,j+1)+U(i,j+1)) - dt/(2*dx)*diff(Fx(:,j+1));
            % Ux(i,j) = (U(i+1,j+1)+U(i,j+1))/2 -  ...
                     % dt/(2*dx)*((U(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                                % (U(i,j+1).^2./H(i,j+1) + g/2*H(i,j+1).^2));
            % y momentum
            Fx = U .* V ./ H;
            Vx = 1/2*(V(i+1,j+1)+V(i,j+1)) - dt/(2*dx)*diff(Fx(:,j+1));
            % Vx(i,j) = (V(i+1,j+1)+V(i,j+1))/2 - ...
                     % dt/(2*dx)*((U(i+1,j+1).*V(i+1,j+1)./H(i+1,j+1)) - ...
                                % (U(i,j+1).*V(i,j+1)./H(i,j+1)));
            % y direction
            i = 1:n;
            j = 1:n+1;
            % height
            Hy = (H(i+1,j+1)+H(i+1,j))/2 - dt/(2*dy)*(V(i+1,j+1)-V(i+1,j));
            % x momentum
            Fy = U.*V./H;
            Uy = 1/2*(U(i+1,j+1)+U(i+1,j)) - dt/(2*dy)*diff(Fy(i+1,:)')';
            % Uy(i,j) = (U(i+1,j+1)+U(i+1,j))/2 - ...
                     % dt/(2*dy)*((V(i+1,j+1).*U(i+1,j+1)./H(i+1,j+1)) - ...
                                % (V(i+1,j).*U(i+1,j)./H(i+1,j)));
            % y momentum
            Fy = (V.^2)./H + g/2*H.^2;
            Vy = 1/2*(V(i+1,j+1)+V(i+1,j)) - dt/(2*dy)*diff(Fy(i+1,:)')';
            % Vy(i,j) = (V(i+1,j+1)+V(i+1,j))/2 - ...
                     % dt/(2*dy)*((V(i+1,j+1).^2./H(i+1,j+1) + g/2*H(i+1,j+1).^2) - ...
                                % (V(i+1,j).^2./H(i+1,j) + g/2*H(i+1,j).^2));
            
            % Second half step
            i = 2:n+1;
            j = 2:n+1;
            % height
            Hnew(i,j) = H(i,j) - (dt/dx)*(Ux(i,j-1)-Ux(i-1,j-1)) - ...
                             (dt/dy)*(Vy(i-1,j)-Vy(i-1,j-1));
            % x momentum
            Fx = (Ux.^2)./Hx + g/2*Hx.^2;
            Fy = Uy.*Vy./Hy;
            Unew(i,j) = U(i,j) - (dt/dx)*diff(Fx) - (dt/dy)*diff(Fy')';
            % Unew(i,j) = U(i,j) - (dt/dx)*((Ux(i,j-1).^2./Hx(i,j-1) + g/2*Hx(i,j-1).^2) - ...
                             % (Ux(i-1,j-1).^2./Hx(i-1,j-1) + g/2*Hx(i-1,j-1).^2)) ...
                           % - (dt/dy)*((Vy(i-1,j).*Uy(i-1,j)./Hy(i-1,j)) - ...
                             % (Vy(i-1,j-1).*Uy(i-1,j-1)./Hy(i-1,j-1)));
            % y momentum
            Fx = Ux.*Vx./Hx;
            Fy = (Vy.^2)./Hy + g/2*Hy.^2;
            Vnew(i,j) = V(i,j) - (dt/dx)*diff(Fx) - (dt/dy)*diff(Fy')';
            % Vnew(i,j) = V(i,j) - (dt/dx)*((Ux(i,j-1).*Vx(i,j-1)./Hx(i,j-1)) - ...
                             % (Ux(i-1,j-1).*Vx(i-1,j-1)./Hx(i-1,j-1))) ...
                           % - (dt/dy)*((Vy(i-1,j).^2./Hy(i-1,j) + g/2*Hy(i-1,j).^2) - ...
                             % (Vy(i-1,j-1).^2./Hy(i-1,j-1) + g/2*Hy(i-1,j-1).^2));
        end % LW
        
        FUx = @(h,u) (u.^2)./h + g/2*h.^2;
        FVx = @(h,u,v) u.*v./h;
        FUy = @(h,u,v) u.*v./h;
        FVy = @(h,v) (v.^2)./h + g/2*h.^2;
        
        if (0)
            % FORCE scheme
            % first half step
            % x
            Hx = 1/2*(H(Is+1,J)+H(Is,J)) - dt/(2*dx) * diff(U(:,J));
            Fx = (U.^2)./H + g/2*H.^2;
            Ux = 1/2*(U(Is+1,J)+U(Is,J)) - dt/(2*dx) * diff(Fx(:,J));
            Fx = U.*V./H;
            Vx = 1/2*(V(Is+1,J)+V(Is,J)) - dt/(2*dx) * diff(Fx(:,J));
            % y
            Hy = 1/2*(H(I,Js+1)+H(I,Js)) - dt/(2*dy) * diff(V(I,:)')';
            Fy = U.*V./H;
            Uy = 1/2*(U(I,Js+1)+U(I,Js)) - dt/(2*dy) * diff(Fy(I,:)')';
            Fy = (V.^2)./H + g/2*H.^2;
            Vy = 1/2*(V(I,Js+1)+V(I,Js)) - dt/(2*dy) * diff(Fy(I,:)')';
            
            % second half step
            Hnew(I,J) = 1/4*(Hx(I-1,:)+Hx(I,:)+Hy(:,J-1)+Hy(:,J)) ...
            - dt/(2*dx)*diff(Ux) - dt/(2*dy)*diff(Uy')';
            %
            Fx = (Ux.^2)./Hx + g/2*Hx.^2;
            Fy = Uy.*Vy./Hy;
            % size(Ux)
            % size(Uy)
            % size(diff(Fx))
            % size(Fy
            % size(diff(Fy')')
            Unew(I,J) = 1/4*(Ux(I-1,:)+Ux(I,:)+Uy(:,J-1)+Uy(:,J)) ...
            - dt/(2*dx)*diff(Fx) - dt/(2*dy)*diff(Fy')';
            %
            Fx = Ux.*Vx./Hx;
            Fy = (Vy.^2)./Hy + g/2*Hy.^2;
            Vnew(I,J) = 1/4*(Vx(I-1,:)+Vx(I,:)+Vy(:,J-1)+Vy(:,J)) ...
            - dt/(2*dx)*diff(Fx) - dt/(2*dy)*diff(Fy')';
        end % FORCE
        
        H(I,J) = Hnew(I,J);
        U(I,J) = Unew(I,J);
        V(I,J) = Vnew(I,J);


        % Update plot
        if mod(nstep,nplotstep) == 0
          C = abs(U(I,J)) + abs(V(I,J));  % Color shows momemtum
          t = nstep*dt;
          tv = norm(C,'fro');
          set(surfplot,'zdata',H(I,J)+B(I,J),'cdata',C);
          set(surfplot,'EdgeColor','none');
          set(bottomplot,'zdata',B(I,J),'cdata',abs(B(I,J)));
          % set(bottomplot,'EdgeColor','none');
          set(top,'string',sprintf('t = %6.2f,  tv = %6.2f',t,tv))
          
          drawnow
        end
        
        if all(all(isnan(H))), break, end  % Unstable, restart
    end % simulation loop
end % main loop
close(gcf)
end
% ------------------------------------

function D = droplet(height,width)
% DROPLET  2D Gaussian
% D = droplet(height,width)
   [x,y] = ndgrid(-1:(2/(width-1)):1);
   D = height*exp(-5*(x.^2+y.^2));
end
% ------------------------------------

function [surfplot,top,restart,quit,bottomplot] = initgraphics(n);
% INITGRAPHICS  Initialize graphics for waterwave.
% [surfplot,top,restart,quit] = initgraphics(n)
% returns handles to a surface plot, its title, and two uicontrol toggles.
    
    clf
    shg
    set(gcf,'numbertitle','off','name','Shallow_water')
    x = (0:n-1)/(n-1);
    surfplot = surf(x,x,ones(n,n),zeros(n,n));
    set(surfplot,'EdgeColor','none');
    grid off
    axis([0 1 0 1 -1 3])
    caxis([-1 1])
    shading faceted
    c = (1:64)'/64;
    % cyan = [0*c c c];
    % colormap(cyan);
    % blue = [0*c, 0*c, c];
    % colormap(blue);
    
    if (1)
        hold on;
        bottomplot = surf(x,x,ones(n,n),zeros(n,n));
        grid off
        axis([0 1 0 1 -1 3])
        % caxis([-1 1])
        shading faceted
        % colormap(jet)
        hold off;
    end
    
    top = title('xxx');
    restart = uicontrol('position',[20 20 80 20],'style','toggle','string','restart');
    quit = uicontrol('position',[120 20 80 20],'style','toggle','string','close');
end


function S = add(A)
S = A(1:end-1,:) + A(2:end,:);
% if (d == 1) 
    % S = A(1:end-1,:) + A(2:end,:);
% elseif (d == 2)
    % S = A(:,1:end-1) + A(:,2:end);
% end
end