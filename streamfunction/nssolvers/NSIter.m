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
	qdbcE2 = filtering.q.outer2.dbc;
	fi = fields(qdbcfullE);
	for i = 1:numel(fi)
		[qdbcEinE2.(fi{i}),qdbcfullEinE2.(fi{i})] = placeInE(qdbcfullE.(fi{i}),filtering.q.outer2.valind,nx+1,ny+1);
	end
	[~,qbciofullEinE2] = placeInE(filtering.q.outer1.bciofull,filtering.q.outer2.valind,nx+1,ny+1);
	
	%biharmonic2(nx,ny,hx,hy,bcxd,bcyd,bcxn,bcyn,posneg)
	%bih = biharmonic2(nx+3,ny+3,hx,hy,qdbcfullE2.w|qdbcfullE2.e,qdbcfullE2.n|qdbcfullE2.s|qdbcfullEinE2.n,(qdbcfullEinE2.w|qdbcfullEinE2.e)&~qdbcfullEinE2.n,qdbcfullEinE2.s,1);
	bih = biharmonic2(nx+3,ny+3,hx,hy,...
		qdbcfullE2.w|qdbcfullE2.e,...
		qdbcfullE2.n|qdbcfullE2.s|(qdbcfullEinE2.n),...
		qdbcfullEinE2.s&~(qdbcfullEinE2.w|qdbcfullEinE2.e),...
		(qdbcfullEinE2.w|qdbcfullEinE2.e)&~qdbcfullEinE2.n,...
		-1);
	%bih = biharmonic2(nx+3,ny+3,hx,hy,[],[],qdbcfullE2.w|qdbcfullE2.e|qdbcfullEinE2.w|qdbcfullEinE2.e&~filtering.q.outer2.bciofull&~qbciofullEinE2,...
	%				qdbcfullE2.s|qdbcfullEinE2.s,1);
	%bih = biharmonic2(nx+3,ny+3,hx,hy,qdbcfullE2.w|qdbcfullE2.e,qdbcfullE2.s|qdbcfullE2.n,[],[],1);
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
	%perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	Dx = sptoeplitz([0 -1],[0 1],nx+3)./(2*hx);
	dx = kron(speye(ny+3),Dx);
	dx = filterMat*dx*filterMat';

	Dy = sptoeplitz([0 -1],[0 1],ny+3)./(2*hy);
	dy = kron(Dy,speye(nx+3));
	dy = filterMat*dy*filterMat';
	
	R = par.Re;
	psi1 = rhs.q.outer2;
	I = speye(size(bih));
	
	rhsbc = rhs.q.outer2;
	%rhsbc(qOnEinE2) = rhs.q.outer1(filtering.q.outer1.on);
	fprintf('\n')
	lents = 0;
	Rnq = RN(psi1,R,qdbcEinE2,qdbcE2,filtering.q.outer2.on|qOnEinE2,rhsbc,dx,dy,Lq,grids,filtering);
	psi1 = bih\Rnq;
	
	for k = 1:par.steps
		%rhs = rhsbc;
		Rnq = RN(psi1,R,qdbcEinE2,qdbcE2,filtering.q.outer2.on|qOnEinE2,rhsbc,dx,dy,Lq,grids,filtering);
		%ivec = grids.q.outer2.yinit(randperm(numel(grids.q.outer2.yinit)));
		psii = bih\Rnq;
		
% 		for i = 1:numel(ivec)
% 			indices = (grids.q.outer2.ymesh==ivec(i));
% 			
% 			M = I;
% 			M(indices,:) = bih(indices,:);
% 			rhs(indices) = Rnq(indices);
% 			psii = M\rhs;
% 			rhs(indices) = (1 - par.omega).*psi1(indices) + par.omega.*psii(indices);
% 			
% 		end

		psi1 = (1 - par.omega).*psi1 + par.omega.*psii;
		
		if((k==1 || mod(k,par.plotoniter)==0)&&~par.noplot)
			res = InPost(psi1,grids,filtering,par);
		end
		
		timestr = sprintf('Step: %d/%d',k,par.steps);
		fprintf([repmat('\b',1,lents) timestr]);
		lents = numel(timestr);
		
	end
	
end

function nq = RN(q,R,qdbcEinE2,qdbcE2,bc,rhs,dx,dy,Lq,grids,filtering)
	%q(bc) = rhs(bc);
	omg = Lq*q;
	u = dy*q;
	v = dx*q;
	u(qdbcEinE2.w|qdbcEinE2.e|qdbcE2.w|qdbcE2.e) = 0;
	v(qdbcEinE2.s|qdbcEinE2.n|qdbcE2.s|qdbcE2.n) = 0;
	omg(bc) = 0;
	
	nq = R*(spdiag(v)*dy*omg-spdiag(u)*dx*omg);
	nq(bc) = rhs(bc);
end

function [inE,fullinE] = placeInE(full,valind,nx,ny)
	fullMat = reshape(full,[nx,ny])';
	fullMatInE = [zeros([ny,1]),fullMat,zeros([ny,1])];
	fullMatInE = logical([zeros([1,nx+2]);fullMatInE;zeros([1,nx+2])]);
	fullinE = reshape(fullMatInE',[],1);
	inE = fullinE(valind);
end

