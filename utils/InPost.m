function [figs,mat,vec] = InPost(qmesh,bc,grids,filtering,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	qmeshfull = filtering.filterMat'*qmesh;
	
	nx = grids.nxp1;
	ny = grids.nyp1;
	
	bcxdfull = bc{1}{2}{1};
	bcydfull = bc{1}{2}{2};
	
	%TODO change so that derivatives are cool at boundaries
	if(par.ghostpoints)

		if(par.filter)
			
			%again we're setting a hard limit of order 3
			switch par.order
				case 1
					%do nothing
				case 2
					[~,~,~,bcxdfull] = closure(grids,filtering,'inner',filtering.gp{1},bcxdfull);
					[~,~,~,bcydfull] = closure(grids,filtering,'inner',filtering.gp{1},bcydfull);
					[~,~,~,qmeshfull] = closure(grids,filtering,'inner',filtering.gp{1},qmeshfull);
					[~,grids,filtering] = closure(grids,filtering,'inner',filtering.gp{1},filtering.gp{1});
				case 3
					[~,~,~,bcxdfull] = closure(grids,filtering,'inner',filtering.gp{2},bcxdfull);
					[~,~,~,bcydfull] = closure(grids,filtering,'inner',filtering.gp{2},bcydfull);
					[~,~,~,qmeshfull] = closure(grids,filtering,'inner',filtering.gp{2},qmeshfull);
					[~,newgrids,newfiltering,gp] = closure(grids,filtering,'inner',filtering.gp{2},filtering.gp{1});

				% 	[~,~,~,umeshfull] = closure(newgrids,newfiltering,h,'inner',gp,qmeshfull);
				% 	[~,~,~,vmeshfull] = closure(newgrids,newfiltering,h,'inner',gp,qmeshfull);
				% 	[~,~,~,qmeshfull] = closure(newgrids,newfiltering,h,'inner',gp,qmeshfull);
				% 	[~,newgrids,newfiltering,gp] = closure(newgrids,newfiltering,h,'inner',gp,gp);

					[~,~,~,bcxdfull] = closure(newgrids,newfiltering,'inner',gp,bcxdfull);
					[~,~,~,bcydfull] = closure(newgrids,newfiltering,'inner',gp,bcydfull);
					[~,~,~,qmeshfull] = closure(newgrids,newfiltering,'inner',gp,qmeshfull);
					[~,grids,filtering] = closure(newgrids,newfiltering,'inner',gp,gp);
				otherwise
					ME = MException('closure:invalidParameterException','Invalid value for par.order');
					throw(ME)
			end
			
		end
		
% 		%TODO figure out how to get back our q at the right size
		filterMat = filtering.filterMat;
		on = filtering.on;
		
		bcxd = logical(filterMat*(1*bcxdfull));
		bcyd = logical(filterMat*(1*bcydfull));
		xmesh = grids.inner.xmesh;
		ymesh = grids.inner.ymesh;
		qmesh = filterMat*qmeshfull;
		nx = grids.nxp1;
		ny = grids.nyp1;
		h = grids.h;
		
		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		dx = spdiag(~(bcxd|bcyd))*dx;
		
		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		dy = spdiag(~(bcxd|bcyd)|(xmesh==min(xmesh)&(ymesh>-.5+h/2&ymesh<0.5-h/2))|(xmesh==max(xmesh)&(ymesh>-1.5+h/2&ymesh<1.5-h/2))...
					     |(xmesh==min(xmesh(on))&(ymesh>-.5+h/2&ymesh<0.5-h/2))|(xmesh==max(xmesh(on))&(ymesh>-1.5+h/2&ymesh<1.5-h/2)))*dy;

		
	else
		filterMat = filtering.filterMat;
		valind = filtering{2};
		on = filtering.on;
		onfull = filtering.onfull;
		
		%Switch to first order on the boundary
		dbc = boundarysides(grids,filtering);
		bcw = dbc.w;
		bce = dbc.e;
		bcs = dbc.s;
		bcn = dbc.n;
		bcc = dbc.c;

		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		
		if(par.zeroout)
			dx = spdiag(~(bcw|bce|bcc))*dx;
		else
			dx = spdiag(~bcw)*dx + 1/h*(-spdiag(bcw) + spdiag(bcw(1:end-1),1));
			dx = spdiag(~bce)*dx + 1/h*(-spdiag(bce(2:end),-1) + spdiag(bce));
		end

		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		
		if(par.zeroout)
			dy = spdiag(~(bcs|bcn|bcc))*dy;
		else
			dy = spdiag(~bcs)*dy + 1/h*(-spdiag(bcs) + spdiag(bcs(1:end-nx),nx));
			dy = spdiag(~bcn)*dy + 1/h*(-spdiag(bcn(nx+1:end),-nx) + spdiag(bcn));
		end
	end
	
	umesh = dy*qmesh;
	vmesh = -dx*qmesh;
	
	umeshfull = filterMat'*umesh;
	Umesh = reshape(umeshfull,[nx,ny])';
	
	vmeshfull = filterMat'*vmesh;
	Vmesh = reshape(vmeshfull,[nx,ny])';
	
	qmeshfull = filterMat'*qmesh;
	Qmesh = reshape(qmeshfull,[nx,ny])';
	
	% use matrices rather than cell arrays so they throw a dimension error if we have a bug mismatched
	mat = cat(3,grids.inner.Xmesh,grids.inner.Ymesh,Umesh,Vmesh,Qmesh);
	vec = cat(2,grids.inner.xmesh,grids.inner.ymesh,umesh,vmesh,qmesh);
	
	if(par.plot)
		if(exist('figs','var'))
			figs = Plot(mat,vec,par,figs);
		else
			figs = Plot(mat,vec,par);
		end
	end
	
end

