
classdef FESpace < handle
    
    properties
        
        mesh;
        
        type;
        vdim;
        
        nDofs;
        nSize;
        
        nNode;
        nElem;
        
        nodes;
        elems;
        
    end % properties
    
    methods
        
        function [ ids ] = getElemNodeIds(fes, ielem)
            ids = fes.elems(ielem,:);
            return
        end
        
        function [ coords ] = getElemNodeCoords(fes, ielem)
            coords = fes.nodes(fes.elems(ielem,:),:);
            return
        end
        
        function [ idofs ] = getElemDofs(fes, ielem)
            nv = fes.vdim;
            ids = fes.elems(ielem,:);
            
            idofs = repmat(ids, nv,1);
            for j = 1:nv
                idofs(j,:) = nv*(ids-1) + j;
            end
            
            idofs = idofs(:);
            
            return
        end
        
        % TODO
        function [ng,xg,wg] = getElemGauss(fes, ielem)
            
            switch fes.type
            case 'Q1'
                nn = 4; % 4 nodes
                ng = 4; % 4 gauss
                [xg,wg] = GaussPoint(nn,ng);
            otherwise
                error('Unsupported gauss for %s', fes.type);
            end
            
            return
        end
        
        function [c] = getCoordSys(fes)
            c = fes.mesh.coord;
            return
        end
        
        
    end % methods
    
    
    methods(Static)
        
    end % methods(Static)
    
    
end





