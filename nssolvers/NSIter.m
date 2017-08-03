function [grids,filtering,res,par,figs] = NSIter(par,grids,filtering,rhs,figs)
	
	filterMat = filtering.outer.filterMat;
	nx = par.nx;
	ny = par.ny;
	h = grids.h;
	
	qOnfullMat = reshape(filtering.inner.onfull,[nx-1,ny-1])';
	qOnfullMatE = [zeros([ny-1,1]),qOnfullMat,zeros([ny-1,1])];
	qOnfullMatE = logical([zeros([1,nx+1]);qOnfullMatE;zeros([1,nx+1])]);
	qOnfullE = reshape(qOnfullMatE',[],1);
	qOnE = logical(filterMat*(1*qOnfullE));
	
	%(nx,ny,h,a11x,a11y,a11io,order,bcwe,bcsn,bcio,posneg)
	bih = biharmonic2(nx+1,ny+1,h,filtering.outer.dbcfull.w|filtering.outer.dbcfull.e,filtering.outer.dbcfull.s|filtering.outer.dbcfull.n,[],[],1);
	bih = filterMat*bih*filterMat';
	
	lap = laplacian2(nx+1,ny+1,h,1,1,1,1,[],[],[],-1);
	lap = filterMat*lap*filterMat';
	
	Dx = sptoeplitz([0 -1],[0 1],nx+1)./(2*h);
	dx = kron(speye(ny+1),Dx);
	dx = filterMat*dx*filterMat';

	Dy = sptoeplitz([0 -1],[0 1],ny+1)./(2*h);
	dy = kron(Dy,speye(nx+1));
	dy = filterMat*dy*filterMat';
	
	qmesh = grids.outer.xmesh*0;
	
	R = 1/par.Re;
	
	rhsbc = rhs.outer;
	%rhsbc(qOnE) = rhs.inner(filtering.inner.on);

	for k = 1:par.timesteps
		Rnq = N(qmesh,R,filtering.outer.on,rhsbc);
		%add a little noise
		%Rnq = Rnq + 0.1*~inout.*rand(numel(inout),1);
		
		qnew = qmesh*0;
		ivec = grids.outer.yinit(randperm(numel(grids.outer.yinit)));
		for i = 1:ny+1
			ycurr = ivec(i);
			indices = (grids.outer.ymesh==ycurr);
			RnqSlice = Rnq(indices);
			bihSlice = (bih(indices,:))';
			bihSlice = (bihSlice(indices,:))';
			
			qnew(indices) = bihSlice\RnqSlice;
			
		end
		
		qmesh = (1-par.omega)*qmesh + par.omega*qnew;

		if(exist('figs','var'))
			[res,figs] = InPost(qmesh,grids,filtering,par,figs);
		else
			[res,figs] = InPost(qmesh,grids,filtering,par);
		end
		
	end
	
	
	function nq = N(q,R,bc,rhs)
		nq = R*spdiag(dy*q)*dx*lap*q-spdiag(dx*q)*dy*lap*q;
		nq(bc) = rhs(bc);
	end
end

