
classdef FEBilinearForm < handle
    properties
        fes;
        mat;
    end
    
    methods
        
        function [febilin] = FEBilinearForm(fespace)
            febilin.fes = fespace;
            febilin.mat = sparse(fespace.nSize,fespace.nSize);
            return
        end
        
        function [] = assembleDomainLaplacian(febilin)
            % the FE space
            fes = febilin.fes;
            
            for ielem = 1:fes.nElem
                iNodePos = fes.getElemNodeCoords(ielem);
                iDofs = fes.getElemDofs(ielem);
                nn = length(iDofs);
                
                % element matrix
                Ke = zeros(nn,nn);
                
                % load Gauss points
                [ngauss,gpoints,weights] = fes.getElemGauss(ielem);
                
                for kGauss = 1:ngauss
                    localCoord = gpoints(kGauss,:);
                    wgt = weights(kGauss);
                    
                    % shape function in local world
                    [N,dN] = LagrangeBasis(fes.type, localCoord);
                    
                    % Jacobian matrix
                    [Jmat,invJ,detJ] = JacobianMatrix(dN, iNodePos);
                    
                    % shape function in real world
                    dN = (Jmat \ dN')';
                    dNdX = dN(:,1);
                    dNdY = dN(:,2);
                    
                    % 
                    if strcmp(fes.getCoordSys(), 'RZ')
                        xg = dot(N,iNodePos(:,1));
                        yg = dot(N,iNodePos(:,2));
                        thick = xg;
                    else
                        thick = 1.0;
                    end
                    
                    % element stiffness
                    Ke = Ke + (dN*dN') * detJ*wgt * thick;
                end
                
                % assemble
                febilin.mat(iDofs,iDofs) = febilin.mat(iDofs,iDofs) + Ke;
            end
            
            return
        end
    end
    
    
end






