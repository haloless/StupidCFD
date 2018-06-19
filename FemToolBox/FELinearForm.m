
classdef FELinearForm < handle
    
    properties
        fes;
        vec;
    end %
    
    methods
        
        function [felin] = FELinearForm(fespace)
            felin.fes = fespace;
            felin.vec = zeros(fespace.nSize, 1);
            return
        end
        
        function [] = setZero(felin)
            felin.vec(:) = 0;
            return
        end
        
        function [] = assembleDomain(felin, fun)
            % the FE space
            fes = felin.fes;
            
            % result vector
            fvec = zeros(size(felin.vec));
            
            for ielem = 1:fes.nElem
                iNodePos = fes.getElemNodeCoords(ielem);
                iDofs = fes.getElemDofs(ielem);
                nn = length(iDofs);
                
                % element vector
                fe = zeros(nn,1);
                
                % load Gauss points
                [ngauss,gpoints,weights] = fes.getElemGauss(ielem);
                
                for kGauss = 1:ngauss
                    localCoord = gpoints(kGauss,:);
                    % xi = gpoints(kGauss,1);
                    % eta = gpoints(kGauss,2);
                    wgt = weights(kGauss);
                    
                    % shape function in local world
                    [N,dN] = LagrangeBasis(fes.type, localCoord);
                    
                    % Jacobian matrix
                    [Jmat,invJ,detJ] = JacobianMatrix(dN, iNodePos);
                    
                    % shape function in real world
                    % dN = (Jmat \ dN')';
                    % dNdX = dN(:,1);
                    % dNdY = dN(:,2);
                    
                    % element stiffness
                    % Ke = Ke + (dN*dN') * detJ*wgt;
                    
                    % external source
                    % global coordinate TODO transform class
                    xg = dot(N,iNodePos(:,1));
                    yg = dot(N,iNodePos(:,2));
                    
                    % TODO 2D/3D
                    fg = fun(xg,yg);
                    fg = reshape(fg(:) * N.', [],1);
                    
                    % 
                    if strcmp(fes.getCoordSys(), 'RZ')
                        thick = xg;
                    else
                        thick = 1.0;
                    end
                    
                    fe = fe + fg * detJ*wgt * thick;
                end
                
                % assemble
                % Kmat(iDofs,iDofs) = Kmat(iDofs,iDofs) + Ke;
                
                fvec(iDofs) = fvec(iDofs) + fe;
            end
            
            % add to linear form
            felin.vec = felin.vec + fvec;
            
            return
        end
        
        
        
    end
    
end


