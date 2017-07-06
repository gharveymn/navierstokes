function [figs,mat,vec] = NSIter(par,figs,mat,vec,grids,filtering,bc,rhs)
	
	qmesh = vec(:,5);
	
	%adjust these
	R = par.Re;
	omega = 0.5;
	
	filterMat = filtering.filterMat;
	nx = grids.nxp1;
	ny = grids.nyp1;
	h = grids.h;
	xinit = grids.inner.xinit;
	yinit = grids.inner.yinit;
	xmesh = grids.inner.xmesh;
	ymesh = grids.inner.ymesh;
	
	bih = biharmonic2(nx,ny,h,bc{1}{2}{1},bc{1}{2}{2});
	bih = filterMat*bih*filterMat';
	
	lap = laplacian2(nx,ny,h);
	lap = filterMat*lap*filterMat';
	
	Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
	dx = kron(speye(ny),Dx);
	dx = filterMat*dx*filterMat';

	Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
	dy = kron(Dy,speye(nx));
	dy = filterMat*dy*filterMat';

	for k = 1:par.timesteps
		Rnq = R*N(qmesh);
		inout = (bc{1}{1}{1}|bc{1}{1}{2});
		Rnq = spdiag(~inout)*Rnq + spdiag(inout)*rhs;
		%add a little noise
		Rnq = Rnq + 0.1*~inout.*rand(numel(inout),1);
		
		qnew = qmesh*0;
% 		for i = 1:ny
% 			ycurr = yinit(i);
% 			indices = (ymesh==ycurr);
% 			RnqSlice = Rnq(indices);
% 			bihSlice = (bih(indices,:))';
% 			bihSlice = (bihSlice(indices,:))';
% 			
% 			qnew(indices) = bihSlice\RnqSlice;
% 			
% 		end

		qnew = bih\Rnq;
		
		[figs,mat,vec] = InPost(qnew,bc,grids,filtering,par,figs);
		
		qmesh = (1-omega)*qnew + omega*qmesh;
	end
	
	
	function nq = N(q)
		nq = spdiag(dy*q)*dx*lap*q-spdiag(dx*q)*dy*lap*q;
	end
end

