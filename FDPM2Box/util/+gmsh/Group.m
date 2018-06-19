%gmsh.Group create a physical group
% 

classdef Group < handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    ndim;
    name;
    conn;
    
    nodeIds;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    
    function [obj] = Group()
        obj.ndim = 0;
        obj.name = '';
        obj.conn = [];
        
        obj.nodeIds = [];
    end
    
    function [ids] = nodes(obj)
        nmax = max(obj.conn(:));
        ids = zeros(nmax, 1);
        for i = 1:size(obj.conn,2)
            ids(obj.conn(:,i)) = 1;
        end
        ids = find(ids == 1);
    end
    
    
end




end






