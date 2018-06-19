%Show eigen vectors of stiffness matrix: INPUT Ktan

if 1
    % check stiffness K matrix eigen structure
    
    figure;
    
    % we copy the K matrix and set Dirichlet BC for boundaries
    % Ktan = Ktan0;
    
    if 1
        % Dirichlet BC
        Ktan(sub2ind([numDofs,numDofs],dirBCDofs,dirBCDofs)) = 1.0e10;
    end
    
    neig = 20;
    % [v,d] = eigs(Ktan, neig);
    [v,d] = eigs(Ktan, neig, 'sm'); % from min eigenvalue
    
    for ieig = 1:neig
        vv = v(:,ieig);
        % surf(xx,yy, reshape(vv,nlen,nlen), 'FaceColor','interp');
        trisurf(tri,nodeX,nodeY, vv, 'FaceColor','interp');
        title(['eig(',int2str(ieig),')=',num2str(d(ieig,ieig))]);
        pause;
    end
end
