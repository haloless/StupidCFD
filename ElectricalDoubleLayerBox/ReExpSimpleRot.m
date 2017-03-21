
%% original Cartesian base
ex = [ 1; 0; 0 ];
ey = [ 0; 1; 0 ];
ez = [ 0; 0; 1 ];

E = [ ex, ey, ez ];

%% the rotated coord

% use center connection as new z-axis
ezhat = pab ./ norm(pab);

exhat = cross(ezhat,ez);
if norm(exhat) == 0
    exhat = ex;
end
exhat = exhat ./ norm(exhat);

eyhat = cross(ezhat,exhat);
eyhat = eyhat ./ norm(eyhat);

Ehat = [ exhat, eyhat, ezhat ];

%% transform matrix for xhat = Q.x
Qmat = Ehat' * E;

% original z-axis in new coord
ezprime = Qmat * ez;
[~,thetaprime,phiprime] = sh_cart2sph(ezprime(1),ezprime(2),ezprime(3));

% new z-axis in orignal coord
[~,anggamma,angchi] = sh_cart2sph(ezhat(1),ezhat(2),ezhat(3));

if 0
    figure;
    plot3([pa(1),pb(1)],[pa(2),pb(2)],[pa(3),pb(3)],'x-')
    % axis equal;
    hold on;
    quiver3(pa(1),pa(2),pa(3), exhat(1),exhat(2),exhat(3));
    quiver3(pa(1),pa(2),pa(3), eyhat(1),eyhat(2),eyhat(3));
    quiver3(pa(1),pa(2),pa(3), ezhat(1),ezhat(2),ezhat(3));
    % quiver3(pb(1),pb(2),pb(3), Ehat(1,:),Ehat(2,:),Ehat(3,:));
    hold off;
    axis equal;
end
