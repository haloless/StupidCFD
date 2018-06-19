%gmsh.PhysicalGroups

classdef PhysicalGroups < handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
properties
    groups;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
methods
    
    function [obj] = PhysicalGroups()
        obj.groups = gmsh.Group();
    end
    
    
    function [g] = subsref(obj, s)
        % s
        % s.type
        % s.subs
        
        ind = s(1).subs{1};
        ind = obj.lookupId(ind);
        g = obj.groups(ind);
        
        if length(s) > 1
            g = subsref(g, s(2:end));
        end
    end
    
    function [obj] = subsasgn(obj, s, b)
        % s
        % stype = s.type
        % ssubs = s.subs
        
        % ind = s.subs{1};
        
        % if ischar(ind)
            
        % else
            % obj.groups(ind) = b;
        % end
        % obj = s;
        
        ind = s(1).subs{1};
        ind = obj.lookupId(ind);
        
        if length(s) > 1
            g = obj.groups(ind);
            g = subsasgn(g, s(2:end), b);
        else
            obj.groups(ind) = b;
        end
    end
    
    function [len] = length(obj)
        len = length(obj.groups);
    end
    
    
    
    function [gid] = getIdByName(obj, name)
        gid = 0;
        for i = 1:length(obj.groups)
            if strcmp(obj.groups(i).name, name)
                gid = i;
                break;
            end
        end
    end
    
    
    
    function [ind] = lookupId(obj, key)
        if ischar(key)
            ind = obj.getIdByName(key);
        elseif isnumeric(key)
            ind = key;
        else
            disp(ind);
            error('Invalid index');
        end
    end
    
end


end




