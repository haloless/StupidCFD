
classdef SpMatAssem < handle
    % SpMatAssem
    % Helper for assemble stiffness matrix as Matlab sparse
    % 
    
    properties
        ii;
        jj;
        ss;
        mm;
        nn;
    end % properties
    
    %
    methods
        
        % construct with matrix size
        function [obj] = SpMatAssem(m,n)
            if nargin <= 0
                error('Must provide matrix size (m,n)');
            elseif nargin == 1
                n = m;
            end
            
            obj.ii = [];
            obj.jj = [];
            obj.ss = [];
            obj.mm = m;
            obj.nn = n;
        end
        
        %
        function [] = SpMatAssemBlockWithDof(obj, Ke,idof,jdof)
            nb = numel(Ke);
            [ib,jb] = ndgrid(idof,jdof);
            obj.ii(end+1:end+nb) = ib(:);
            obj.jj(end+1:end+nb) = jb(:);
            obj.ss(end+1:end+nb) = Ke(:);
        end
        
        %
        function [K] = SpMatCreateAndClear(obj)
            K = sparse(obj.ii,obj.jj,obj.ss, obj.mm,obj.nn);
            
            % clear to save memory and for next use
            obj.ii = [];
            obj.jj = [];
            obj.ss = [];
        end
    end % methods
    
    
end % SpMatAssem








