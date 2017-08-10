function [grids,filtering,res,par] = NSIter(par,grids,filtering,rhs)
	
	if(~exist('par','var'))
		Main();
		return;
	end
	
	filterMat = filtering.q.outer2.filterMat;
	nx = par.nx;
	ny = par.ny;
	h = par.h;
	hx = par.hx;
	hy = par.hy;
	
	qOnfullMatE = reshape(filtering.q.outer1.onfull,[nx+1,ny+1])';
	qOnfullMatE2 = reshape(filtering.q.outer2.onfull,[nx+3,ny+3])';
	qOnfullMatEinE2 = [zeros([ny+1,1]),qOnfullMatE,zeros([ny+1,1])];
	qOnfullMatEinE2 = logical([zeros([1,nx+3]);qOnfullMatEinE2;zeros([1,nx+3])]);
	qOnfullEinE2 = reshape(qOnfullMatEinE2',[],1);
	qOnEinE2 = qOnfullEinE2(filtering.q.outer2.valind);
	
	qdbcfullE = filtering.q.outer1.dbcfull;
	qdbcfullE2 = filtering.q.outer2.dbcfull;
	fi = fields(qdbcfullE);
	for i = 1:numel(fi)
		[~,qdbcfullEinE2.(fi{i})] = placeInE(qdbcfullE.(fi{i}),filtering.q.outer2.valind,nx+1,ny+1);
	end
	
	
	%(nx,ny,h,a11x,a11y,a11io,order,bcwe,bcsn,bcio,posneg)
	bih = biharmonic2(nx+3,ny+3,hx,hy,qdbcfullE2.w|qdbcfullE2.e|qdbcfullEinE2.w|qdbcfullE2.e,qdbcfullE2.s|qdbcfullE2.n|qdbcfullEinE2.s|qdbcfullE2.n,[],[],1);
	bih = filterMat*bih*filterMat';
	
	%stencil boundary conds
	neumannbd = 1;
	dirichletbd = 2;
	dirichletmidbd = 3;
	
	qdbcfull = filtering.q.outer2.dbcfull;
	qbcpar.we.inds = qdbcfull.w|qdbcfull.e;
	qbcpar.sn.inds = qdbcfull.s|qdbcfull.n;
	qbcpar.io.inds = filtering.q.outer2.bciofull;
	
	qbcpar.we.a11.x = dirichletbd;
	qbcpar.we.a11.y = dirichletbd;
	
	qbcpar.sn.a11.x = dirichletbd;
	qbcpar.sn.a11.y = dirichletbd;
	
	qbcpar.io.a11.x = dirichletbd;
	qbcpar.io.a11.y = dirichletbd;
	Lq = laplacian2(nx+3,ny+3,hx,hy,1,-1,qbcpar);
	Lq = filtering.q.outer2.filterMat*Lq*filtering.q.outer2.filterMat';
	perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	Dx = sptoeplitz([0 -1],[0 1],nx+3)./(2*hx);
	dx = kron(speye(ny+3),Dx);
	dx = filterMat*dx*filterMat';

	Dy = sptoeplitz([0 -1],[0 1],ny+3)./(2*hy);
	dy = kron(Dy,speye(nx+3));
	dy = filterMat*dy*filterMat';
	
	qmesh = grids.q.outer2.xmesh*0;
	
	R = 1/par.Re;
	
	rhsbc = rhs.q.outer2;
	rhsbc(qOnEinE2) = rhs.q.outer1(filtering.q.outer1.on);

	
	fprintf('\n')
	lents = 0;
	for k = 1:par.steps
		Rnq = N(qmesh,R,filtering.q.outer2.on|qOnEinE2,rhsbc,dx,dy,Lq);
		%add a little noise
		%Rnq = Rnq + 0.1*~inout.*rand(numel(inout),1);
		
		qnew = qmesh*0;
		ivec = grids.q.outer2.yinit(randperm(numel(grids.q.outer2.yinit)));
% 		for i = 1:numel(ivec)
% 			ycurr = ivec(i);
% 			indices = (grids.q.outer2.ymesh==ycurr);
% 			RnqSlice = Rnq(indices);
% 			bihSlice = (bih(indices,:))';
% 			bihSlice = (bihSlice(indices,:))';
% 			
% 			qnew(indices) = bihSlice\RnqSlice;
% 			
% 		end

		qnew = bih\Rnq;
		
		qmesh = (1-par.omega)*qmesh + par.omega*qnew;
		
		if((k==1 || mod(k,par.plotoniter)==0)&&~par.noplot)
			res = InPost(qmesh,grids,filtering,par);
		end
		
		timestr = sprintf('Step: %d/%d',k,par.steps);
		fprintf([repmat('\b',1,lents) timestr]);
		lents = numel(timestr);
		
	end
	
end

function nq = N(q,R,bc,rhs,dx,dy,Lq)
	q(bc) = rhs(bc);
	nq = R*(spdiag(dy*q)*dx*Lq*q-spdiag(dx*q)*dy*Lq*q);
	nq(bc) = rhs(bc);
end

function [inE,fullinE] = placeInE(full,valind,nx,ny)
	fullMat = reshape(full,[nx,ny])';
	fullMatInE = [zeros([ny,1]),fullMat,zeros([ny,1])];
	fullMatInE = logical([zeros([1,nx+2]);fullMatInE;zeros([1,nx+2])]);
	fullinE = reshape(fullMatInE',[],1);
	inE = fullinE(valind);
end

