function res = InPost(qmesh,grids,filtering,par)
	%INPOST does the post processing of calculation
	
	h = par.h;
	hx = par.hx;
	hy = par.hy;
	
	nx = par.nx;
	ny = par.ny;
	
	filterMat = filtering.q.outer2.filterMat;
	
	%Switch to first order on the boundary
	dbcouter2 = filtering.q.outer2.dbc;
	bcw = dbcouter2.w;
	bce = dbcouter2.e;
	bcs = dbcouter2.s;
	bcn = dbcouter2.n;
	bcc = dbcouter2.c;
	
	Dx = sptoeplitz([0 -1],[0 1],nx+3)./(2*hx);
	dx = kron(speye(ny+3),Dx);
	dx = filterMat*dx*filterMat';
	
	dx = spdiag(~bcw)*dx + 1/hx*(-spdiag(bcw) + spdiag(bcw(1:end-1),1));
	dx = spdiag(~bce)*dx + 1/hx*(-spdiag(bce(2:end),-1) + spdiag(bce));
	
	Dy = sptoeplitz([0 -1],[0 1],ny+3)./(2*hy);
	dy = kron(Dy,speye(nx+3));
	dy = filterMat*dy*filterMat';
	
	dy = spdiag(~bcs)*dy + 1/hy*(-spdiag(bcs) + spdiag(bcs(1:end-nx-4),nx+4));
	dy = spdiag(~bcn)*dy + 1/hy*(-spdiag(bcn(nx+4:end),-nx-3) + spdiag(bcn));
	
	umesh = dy*qmesh;
	vmesh = -dx*qmesh;
	
	umeshfull = filterMat'*umesh;
	res.U = reshape(umeshfull,[nx+3,ny+3])';
	
	vmeshfull = filterMat'*vmesh;
	res.V = reshape(vmeshfull,[nx+3,ny+3])';
	
	qmeshfull = filterMat'*qmesh;
	res.Q = reshape(qmeshfull,[nx+3,ny+3])';
	
	Plot(grids,filtering,res,par);
	
end

