
clear all;

bc_dir = 1;
bc_neu = 2;

% input_file = 'input.dat'; 
% input_file = 'input_ring_M360.dat'; 
input_file = 'input_square_M400q.dat'; 
fid = fopen(input_file);

% problem title
tline = fgetl(fid);
title = tline;
disp(title);

% number
tline = fgetl(fid);
cc = textscan(tline, '%d %d');
n = cc{1};
nfield = cc{2};
% aa = fscanf(

% allocate
% elements
u = zeros(n,1);
x = zeros(2,n);
y = zeros(2,n);
node = zeros(2,n);
bc = zeros(2,n);
dlen = zeros(n,1);
dnorm = zeros(2,n);
% fields
xfield = zeros(2,nfield);

%
% load elements
%

% skip section header
tline = fgetl(fid)
% element end point
for i = 1:n
    tline = fgetl(fid);
    cc = textscan(tline, '%d %f %f');
    y(1,i) = cc{2};
    y(2,i) = cc{3};
end

% skip section header
tline = fgetl(fid)
% element connectivity and BC
for i = 1:n
    tline = fgetl(fid);
    cc = textscan(tline, '%d %d %d %d %f');
    node(1,i) = cc{2};
    node(2,i) = cc{3};
    bc(1,i) = cc{4};
    bc(2,i) = cc{5};
end

if nfield > 0
% skip section header
tline = fgetl(fid)
% field points
for i = 1:nfield
    tline = fgetl(fid);
    cc = textscan(tline, '%d %f %f');
    xfield(1,i) = cc{2};
    xfield(2,i) = cc{3};
end
end

fclose(fid);
clear fid tline cc;

% compute midpoints and normals for elements
for i = 1:n
    j1 = node(1,i);
    j2 = node(2,i);
    x(1,i) = 0.5 * (y(1,j1) + y(1,j2));
    x(2,i) = 0.5 * (y(2,j1) + y(2,j2));
    
    h1 = y(2,j2) - y(2,j1);
    h2 = -y(1,j2) + y(1,j1);
    el = sqrt(h1^2 + h2^2);
    dnorm(1,i) = h1/el;
    dnorm(2,i) = h2/el;
    dlen(i) = el;
end

% plot elements
if (1)
    figure;
    hold on;
    for i = 1:n
        j1 = node(1,i);
        j2 = node(2,i);
        plot([y(1,j1),y(1,j2)],[y(2,j1),y(2,j2)],'.-b');
        text(x(1,i),x(2,i),num2str(bc(:,i)));
    end
    plot(x(1,:),x(2,:),'.r');
    % quiver(x(1,:),x(2,:),dnorm(1,:),dnorm(2,:));
    plot(xfield(1,:),xfield(2,:),'.g');
    hold off;
    axis equal;
end

%
% b-vector
%
bvec = zeros(n,1);
for j = 1:n
    % element point
    al = dlen(j);
    j1 = node(1,j);
    j2 = node(2,j);
    
    for i = 1:n
        % source point
        x11 = y(1,j1) - x(1,i);
        x21 = y(2,j1) - x(2,i);
        x12 = y(1,j2) - x(1,i);
        x22 = y(2,j2) - x(2,i);
        
        r1 = sqrt(x11^2 + x21^2);
        r2 = sqrt(x12^2 + x22^2);
        d = x11*dnorm(1,j) + x21*dnorm(2,j);
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j);
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j);
        
        ds = abs(d);
        theta1 = atan2(t1,ds);
        theta2 = atan2(t2,ds);
        dtheta = theta2 - theta1;
        
        aa = (-dtheta*ds + al + t1*log(r1) - t2*log(r2)) / (pi*2);
        if d < 0
            dtheta = -dtheta;
        end
        bb = -dtheta / (pi*2);
        if i == j
            bb = 0.5;
        end
        
        if bc(1,j) == 1
            bvec(i) = bvec(i) - bb*bc(2,j);
        elseif bc(1,j) == 2
            bvec(i) = bvec(i) + aa*bc(2,j);
        end
    end
end

%
%
%
Amat = zeros(n,n);
for j = 1:n
    al = dlen(j);
    j1 = node(1,j);
    j2 = node(2,j);
    
    for i = 1:n
        x11 = y(1,j1) - x(1,i);
        x21 = y(2,j1) - x(2,i);
        x12 = y(1,j2) - x(1,i);
        x22 = y(2,j2) - x(2,i);
        
        r1 = sqrt(x11^2 + x21^2);
        r2 = sqrt(x12^2 + x22^2);
        d = x11*dnorm(1,j) + x21*dnorm(2,j);
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j);
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j);
        
        ds = abs(d);
        theta1 = atan2(t1,ds);
        theta2 = atan2(t2,ds);
        dtheta = theta2 - theta1;
        
        aa = (-dtheta*ds + al + t1*log(r1) - t2*log(r2)) / (pi*2);
        if d < 0
            dtheta = -dtheta;
        end
        bb = -dtheta / (pi*2);
        
        if i ~= j
            if bc(1,j) == 1
                Amat(i,j) = Amat(i,j) - aa;
            elseif bc(1,j) == 2
                Amat(i,j) = Amat(i,j) + bb;
            end
        else
            if bc(1,j) == 1
                Amat(i,j) = Amat(i,j) - aa;
            elseif bc(1,j) == 2
                Amat(i,j) = Amat(i,j) + 0.5;
            end
        end
    end
end

% solution
sol = Amat \ bvec;
u = sol;

% u contains both value and derivative
% replace all Dirichlet BC
uval = sol;
for i = 1:n
    if bc(1,i) == bc_dir
        uval(i) = bc(2,i);
    end
end



%
% evaluate fields
%
f = zeros(nfield,1);
for j = 1:n
    if bc(1,j) == bc_dir
        f0 = bc(2,j);
        df0 = u(j);
    elseif bc(1,j) == bc_neu
        f0 = u(j);
        df0 = bc(2,j);
    end
    
    al = dlen(j);
    j1 = node(1,j);
    j2 = node(2,j);
    
    for i = 1:nfield
        x11 = y(1,j1) - xfield(1,i);
        x21 = y(2,j1) - xfield(2,i);
        x12 = y(1,j2) - xfield(1,i);
        x22 = y(2,j2) - xfield(2,i);
        
        r1 = sqrt(x11^2 + x21^2);
        r2 = sqrt(x12^2 + x22^2);
        d = x11*dnorm(1,j) + x21*dnorm(2,j);
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j);
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j);
        
        ds = abs(d);
        theta1 = atan2(t1,ds);
        theta2 = atan2(t2,ds);
        dtheta = theta2 - theta1;
        
        aa = (-dtheta*ds + al + t1*log(r1) - t2*log(r2)) / (pi*2);
        if d < 0
            dtheta = -dtheta;
        end
        bb = -dtheta / (pi*2);
        
        f(i) = f(i) + aa*df0 - bb*f0;
    end
end

if (1)
    figure;
    hold on;
    plot3(x(1,:),x(2,:),uval,'x-');
    plot3(xfield(1,:),xfield(2,:),f,'o');
    hold off;
end





