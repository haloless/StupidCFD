
clear all;

% initial grid length
% reflen = 0.2; 
reflen = 0.3; 

% distance function
fd=@(p) dsphere(p,0,0,0,1);

%
[pt0,tri0] = distmeshsurface(fd,@huniform,reflen,1.1*[-1,-1,-1;1,1,1]);

% then generate neighbor list for each element

np = size(pt0,1);
ne = size(tri0,1);


p2e = zeros(np,16);
for ke = 1:ne
    for i = 1:3
        ip = tri0(ke,i);
        p2e(ip,end) = p2e(ip,end)+1;
        if p2e(ip,end)>=size(p2e,2)
            error('p2e overflow')
        end
        p2e(ip,p2e(ip,end)) = ke;
    end
end

nbr0 = zeros(ne,3);
for ke = 1:ne
    for i = 1:3
        i1 = mod(i,3)+1;
        ip = tri0(ke,i);
        ip1 = tri0(ke,i1);
        
        neigh_elem = 0;
        for j = 1:p2e(ip,end)
            je = p2e(ip,j);
            if je==ke; continue; end
            
            match = 0;
            if tri0(je,1)==ip1 && tri0(je,2)==ip; match=1; end
            if tri0(je,2)==ip1 && tri0(je,3)==ip; match=1; end
            if tri0(je,3)==ip1 && tri0(je,1)==ip; match=1; end
            
            if match
                neigh_elem = je;
                break;
            end
        end
        
        if ~neigh_elem
            fprintf('ELEM%d NODE(%d,%d) no match\n',ke,ip,ip1);
            error('Failed to find neighbor element');
        end
        
        nbr0(ke,i) = neigh_elem;
    end
    if mod(ke,50)==0 || ke==ne
        disp(['ke=',int2str(ke)]);
    end
end

if (1)
    cen = zeros(ne,3);
    for ke = 1:ne
        for i = 1:3
            cen(ke,:) = cen(ke,:) + pt0(tri0(ke,i),:);
        end
    end
    cen = cen ./ 3.0;
    
    label = int2str((1:ne)');
    
    figure;
    trimesh(tri0,pt0(:,1),pt0(:,2),pt0(:,3));
    % hidden off;
    axis equal;
    
    hold on
    text(cen(:,1),cen(:,2),cen(:,3), label);
    hold off
end

% dump to file
if (1)
    save('unit_sphere.mat','reflen','pt0','tri0','nbr0');
end

