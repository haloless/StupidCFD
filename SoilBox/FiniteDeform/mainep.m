
clear;

NRitmax = 500000;
tNRits = 0;
gNRits = 0;
NRtol = 1.0e-12;

bm1 = [ 1, 1, 1, 0, 0, 0 ]';

%
[nels,nodes,ndof,coord,etopol,fext0,bc,nen,ngp,lstps,E,v,fc] = setup_cube();
disp(['nelem=',int2str(nels)]);
disp(['node per elem=',int2str(nen)]);
disp(['ngauss=',int2str(ngp)]);

nen3 = nen * 3;

%
[wp,xsi,eta,zet,dNr] = der_shape_func(ngp);

% set element dof
edof = zeros(nels,nen3+1);
dof = [(1:3:ndof-2)', (2:3:ndof-1)', (3:3:ndof)'];
for nel = 1:nels
	% each row for element dof [ielem, node1x,node1y,node1z, ... ]
	edof(nel,1) = nel;
	for node = 1:nen
		edof(nel,3*node-1:3*node+1) = dof(etopol(nel,node),:);
	end
end

bcdof = bc(:,1);
bcval = bc(:,2);


% global tangent stiffness
Ktan0 = zeros(ndof);
% internal force
fint = zeros(ndof,1);
%
rct = zeros(ndof,1);

uvwold = zeros(ndof,1);
uvw = zeros(ndof,1);

% log strain
epsEold = zeros(6,1, nels,ngp);
epsE = epsEold;
% cauchy
sig = zeros(6,1, nels,ngp);
% deform grad
Fold = zeros(3,3,nels,ngp);

%
[ex,ey,ez] = xsplit(edof,coord,dof,nen,nels);

% initial D-matrix
D = vmconst(zeros(6,1), E,v,fc);

% initial deform-grad
for nel = 1:nels
	ke = zeros(nen*3);
	
	% Jacobian dx/deta
	JT = dNr * [ex(nel,:); ey(nel,:); ez(nel,:)]';
	
	for gp = 1:ngp
		Fold(1:3,1:3,nel,gp) = eye(3);
		indx = [3*gp-2; 3*gp-1; 3*gp];
		
		% chain rule, back-solve by Jacobian matrix
		B = formBG(((JT(indx,:)) \ dNr(indx,:)), nen);
		
		% element stiffness
		ke = ke + (B'*D*B*det(JT(indx,:))*wp(gp));
	end
	
	Ktan0 = assem(edof(nel,:),Ktan0,ke);
end

NRtol = max([abs(fext0); abs(Ktan0(bc(:,1),bc(:,1))*bc(:,2))]) * NRtol;
NRtol = max([NRtol, 1.0e-6]);

k = 0;
F = Fold;

% incremental loading steps
for lstp = 1:lstps
	% external loading 
	fext = (lstp/lstps) * fext0;
	
	% out-of-balance force
	oobf = rct + fext - fint;
	
	obfN = 2 * NRtol;
	NRit = 0;
	
	% Newton-Raphson iteration
	while (NRit<NRitmax) && (obfN>NRtol)
		NRit = NRit + 1;
		tNRits = tNRits + 1;
		gNRits = max(gNRits,NRit);
		
		if (NRit == 1)
			[dduvw,dreact] = solveq(Ktan0,oobf,[bcdof, (1/lstps)*bcval]);
		else
			[dduvw,dreact] = solveq(Ktan, oobf,[bcdof, zeros(size(bcval))]);
		end
		
		uvw = uvw + dduvw;
		rct = rct + dreact;
		duvw = uvw - uvwold;
		
		% split nodal displacement into elements
		[eu,ev,ew] = xsplit(edof,[uvw(dof(:,1)),uvw(dof(:,2)),uvw(dof(:,3))],dof,nen,nels);
		[du,dv,dw] = xsplit(edof,[duvw(dof(:,1)),duvw(dof(:,2)),duvw(dof(:,3))],dof,nen,nels);
		
		%
		Ktan = zeros(ndof);
		fint = zeros(ndof,1);
		felem = zeros(nen*3,nels);
		
		for nel = 1:nels
			% Jacobian dx/deta
			JT = dNr * [ex(nel,:)+eu(nel,:); ey(nel,:)+ev(nel,:); ez(nel,:)+ew(nel,:)]';
			
			% element stiffness
			ke = zeros(nen3);
			for gp = 1:ngp
				indx = [3*gp-2; 3*gp-1; 3*gp];
				detJ = det(JT(indx,:));
				dNx = JT(indx,:) \ dNr(indx,:);
				[B,G] = formBG(dNx,nen);
				
				ddF = zeros(3);
				for i = 1:size(dNx,2)
					ddF = ddF + [du(nel,i); dv(nel,i); dw(nel,i)]* dNx(:,i)';
				end
				
				% update deform-grad
				dF = inv(eye(3)-ddF);
				F(:,:,nel,gp) = dF * Fold(:,:,nel,gp);
				
				e = epsEold(:,:,nel,gp);
				eps = [e(1),e(4)/2,e(6)/2; e(4)/2,e(2),e(5)/2; e(6)/2,e(5)/2,e(3)];
				
				% polar decomp
				[V,ePr] = eig(eps);
				ePr = diag(ePr);
				Be = exp(2*ePr);
				Be = Be(1)*V(:,1)*V(:,1)' + Be(2)*V(:,2)*V(:,2)' + Be(3)*V(:,3)*V(:,3)';
				
				% trial b
				BeTr = dF * Be * dF';
				
				[V,BePr] = eig(BeTr);
				BePr = diag(BePr);
				BePr = 0.5 * log(BePr);
				
				% trial e
				etr = BePr(1)*V(:,1)*V(:,1)' + BePr(2)*V(:,2)*V(:,2)' + BePr(3)*V(:,3)*V(:,3)';
				% voigt form
				epsEtr = [etr(1), etr(5), etr(9), 2*etr(2), 2*etr(6), 2*etr(7)]';
				
				% constitutive equation
				% algo tangent, kirchhoff, 
				[D,kirSig,epsE(:,1,nel,gp)] = vmconst(epsEtr,E,v,fc);
				
				% spatial tangent operator [a] and updated cauchy
				[a,sig(:,1,nel,gp)] = formDsig(BeTr, kirSig, D,F(:,:,nel,gp));
				
				ke = ke + (G'*a*G) * detJ*wp(gp);
				felem(:,nel) = felem(:,nel) + B'*sig(:,:,nel,gp)*detJ*wp(gp);
			end
			
			% element force
			eind = edof(nel,2:nen3+1);
			fint(eind) = fint(eind) + felem(:,nel);
			% global stiffness
			Ktan = assem(edof(nel,:),Ktan,ke);
		end
		
		oobf = fext + rct - fint;
		obfN = norm(oobf);
	end
	
	% update state
	uvwold = uvw;
	epsEold = epsE;
	Fold = F;
end





















