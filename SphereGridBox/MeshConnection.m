function [edge,face2edge,edge2face] = MeshConnection(face)

% check face, this can be directly used by Matlab's TRIMESH and TRISURF, etc.
if size(face,2) ~= 3
	error('face must be nface*3');
end

% face number
nface = size(face,1);

%
% build edge
%
nedge = 0;
edge = zeros(nface*3/2, 2);

for i = 1:nface
	i1 = face(i,1);
	i2 = face(i,2);
	i3 = face(i,3);
	
	% i1,i2
	if FindEdge(nedge,edge,i1,i2) == 0
		nedge = nedge + 1;
		if i1 < i2
			edge(nedge,:) = [i1,i2];
		else
			edge(nedge,:) = [i2,i1];
		end
	end
	
	% i2,i3
	if FindEdge(nedge,edge,i2,i3) == 0
		nedge = nedge + 1;
		if i2 < i3
			edge(nedge,:) = [i2,i3];
		else
			edge(nedge,:) = [i3,i2];
		end
	end
	
	% i3,i1
	if FindEdge(nedge,edge,i3,i1) == 0
		nedge = nedge + 1;
		if i3 < i1
			edge(nedge,:) = [i3,i1];
		else
			edge(nedge,:) = [i1,i3];
		end
	end
end

if nedge ~= nface*3/2
	error('edge number not match');
end

% sort the edges by two nodes
% although this is not necessary
edge = sortrows(edge,[1,2]);

%
% match edge and face 
%

% face->edge, edge has direction (sign)
face2edge = zeros(nface,3);
% edge->face, 1st same direction, 2nd inverse direction
edge2face = zeros(nedge,2);

for iface = 1:nface
    i1 = face(iface,1);
    i2 = face(iface,2);
    i3 = face(iface,3);
    
    for jedge = 1:nedge 
        ja = edge(jedge,1);
        jb = edge(jedge,2);
        
        % (i1,i2) 
        if i1==ja && i2==jb
            face2edge(iface,1) = jedge;
            edge2face(jedge,1) = iface;
        elseif i1==jb && i2==ja
            face2edge(iface,1) = -jedge;
            edge2face(jedge,2) = iface;
        end
        
        % i2,i3
        if i2==ja && i3==jb
            face2edge(iface,2) = jedge;
            edge2face(jedge,1) = iface;
        elseif i2==jb && i3==ja
            face2edge(iface,2) = -jedge;
            edge2face(jedge,2) = iface;
        end
        
        % i3,i1
        if i3==ja && i1==jb
            face2edge(iface,3) = jedge;
            edge2face(jedge,1) = iface;
        elseif i3==jb && i1==ja
            face2edge(iface,3) = -jedge;
            edge2face(jedge,2) = iface;
        end
    end
end


return
end

function [jedge] = FindEdge(nedge,edge,ja,jb)
	% ensure a<b
	if ja > jb
		tmp = ja;
		ja = jb;
		jb = tmp;
	end
	
	jedge = 0;
	for j = 1:nedge
		if edge(j,1)==ja && edge(j,2)==jb
			jedge = j;
			break;
		end
	end
	
	return
end



