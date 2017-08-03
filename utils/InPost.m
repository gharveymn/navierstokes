function [res,figs] = InPost(qmesh,grids,filtering,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	qmeshfull = filtering.outer.filterMat'*qmesh;
	
	nx = par.nx;
	ny = par.ny;
	
	%TODO change so that derivatives are cool at boundaries
	if(par.ghostpoints)
		
% 		%TODO figure out how to get back our q at the right size
		filterMat = filtering.outer.filterMat;
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
		filterMat = filtering.outer.filterMat;
		
		%Switch to first order on the boundary
		dbc = boundarysides(grids,filtering,par,'outer',nx);
		bcw = dbc.w;
		bce = dbc.e;
		bcs = dbc.s;
		bcn = dbc.n;
		bcc = dbc.c;

		Dx = sptoeplitz([0 -1],[0 1],nx+1)./(2*h);
		dx = kron(speye(ny+1),Dx);
		dx = filterMat*dx*filterMat';
		
		if(par.zeroout)
			dx = spdiag(~(bcw|bce|bcc))*dx;
		else
			dx = spdiag(~bcw)*dx + 1/h*(-spdiag(bcw) + spdiag(bcw(1:end-1),1));
			dx = spdiag(~bce)*dx + 1/h*(-spdiag(bce(2:end),-1) + spdiag(bce));
		end

		Dy = sptoeplitz([0 -1],[0 1],ny+1)./(2*h);
		dy = kron(Dy,speye(nx+1));
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
	res.U = reshape(umeshfull,[nx+1,ny+1])';
	
	vmeshfull = filterMat'*vmesh;
	res.V = reshape(vmeshfull,[nx+1,ny+1])';
	
	qmeshfull = filterMat'*qmesh;
	res.Q = reshape(qmeshfull,[nx+1,ny+1])';
	
	if(par.plot)
		if(exist('figs','var'))
			figs = Plot(grids,filtering,res,par,figs);
		else
			figs = Plot(grids,filtering,res,par);
		end
	end
	
end

