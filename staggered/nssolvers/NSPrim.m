function [grids,filtering,res,par] = NSPrim(par,grids,filtering,rhs)
	%NSPRIM Strang primitive formulation of NS
	
	if(~exist('par','var'))
		Main();
		return;
	end
	
	nx = par.nx;
	ny = par.ny;
	hx = par.hx;
	hy = par.hy;
	
	U = zeros([ny,nx-1]);
	V = zeros([ny-1,nx]);
	Ue = zeros([ny+2,nx+1]);
	Ve = zeros([ny+1,nx+2]);
	
	if(par.useGPU)
		U = gpuArray(U);
		V = gpuArray(V);
		Ue = gpuArray(Ue);
		Ve = gpuArray(Ve);
	end
	
	%inner
	udbcfull = filtering.u.inner.dbcfull;
	vdbcfull = filtering.v.inner.dbcfull;
	pdbcfull = filtering.p.inner.dbcfull;
	qdbcfull = filtering.q.inner.dbcfull;
	
	%outer
	udbcfullE = filtering.u.outer.dbcfull;
	vdbcfullE = filtering.v.outer.dbcfull;
	pdbcfullE = filtering.p.outer.dbcfull;
	qdbcfullE = filtering.q.outer.dbcfull;
	
	for i=1:numel(par.varnames)
		matlab.lang.makeValidName([par.varnames{i} 'dbcfull']);
		eval([par.varnames{i} 'dbcfull = filtering.' par.varnames{i} '.inner.dbcfull;']);
		eval([par.varnames{i} 'dbcfullE = filtering.' par.varnames{i} '.outer.dbcfull;']);
	end
	
	%gridify
	dbcfields = fields(udbcfull);
	
	for j=1:numel(dbcfields)
		
		Udbcfull.(dbcfields{j}) = gridify(udbcfull.(dbcfields{j}),nx-1,ny);
		Vdbcfull.(dbcfields{j}) = gridify(vdbcfull.(dbcfields{j}),nx,ny-1);
		Pdbcfull.(dbcfields{j}) = gridify(pdbcfull.(dbcfields{j}),nx,ny);
		Qdbcfull.(dbcfields{j}) = gridify(qdbcfull.(dbcfields{j}),nx-1,ny-1);
		
		UdbcfullE.(dbcfields{j}) = gridify(udbcfullE.(dbcfields{j}),nx+1,ny+2);
		VdbcfullE.(dbcfields{j}) = gridify(vdbcfullE.(dbcfields{j}),nx+2,ny+1);
		PdbcfullE.(dbcfields{j}) = gridify(pdbcfullE.(dbcfields{j}),nx+2,ny+2);
		QdbcfullE.(dbcfields{j}) = gridify(qdbcfullE.(dbcfields{j}),nx+1,ny+1);
		
	end
	
	uselEw = UdbcfullE.w&~(UdbcfullE.s|UdbcfullE.n);
	uselEe = UdbcfullE.e&~(UdbcfullE.s|UdbcfullE.n);
	%uselEs = UdbcfullE.s&~(UdbcfullE.w|UdbcfullE.e);
	%uselEn = UdbcfullE.n&~(UdbcfullE.w|UdbcfullE.e);
	
	%vselEw = VdbcfullE.w&~(VdbcfullE.s|VdbcfullE.n);
	%vselEe = VdbcfullE.e&~(VdbcfullE.s|VdbcfullE.n);
	vselEs = VdbcfullE.s&~(VdbcfullE.w|VdbcfullE.e);
	vselEn = VdbcfullE.n&~(VdbcfullE.w|VdbcfullE.e);
	
	%inner boundary cond
	Urhs = reshape(filtering.u.inner.filterMat'*rhs.u.inner,[nx-1,ny])';
	Ubc = Urhs;
	Ubc(Udbcfull.n|Udbcfull.s) = 2*par.dt/(hx^2*par.Re)*Urhs(Udbcfull.n|Udbcfull.s);
	Ubc(Udbcfull.w|Udbcfull.e) = par.dt/(hy^2*par.Re)*Urhs(Udbcfull.w|Udbcfull.e);
	
	Vrhs = reshape(filtering.v.inner.filterMat'*rhs.v.inner,[nx,ny-1])';
	Vbc = Vrhs;
	Vbc(Vdbcfull.n|Vdbcfull.s) = par.dt/(hx^2*par.Re)*Vrhs(Vdbcfull.s|Vdbcfull.n);
	Vbc(Vdbcfull.w|Vdbcfull.e) = 2*par.dt/(hy^2*par.Re)*Vrhs(Vdbcfull.w|Vdbcfull.e);
	
	QUbc = 0*reshape(filtering.q.inner.filterMat'*rhs.q.inner,[nx-1,ny-1])';
	QVbc = QUbc;
	
	%outer boundary cond
	UbcE = reshape(filtering.u.outer.filterMat'*rhs.u.outer,[nx+1,ny+2])';
	VbcE = reshape(filtering.v.outer.filterMat'*rhs.v.outer,[nx+2,ny+1])';
	QUbcE = 0*reshape(filtering.q.outer.filterMat'*rhs.q.outer,[nx+1,ny+1])';
	QVbcE = QUbcE;
	
	%stencil boundary conds
	neumannbd = 1;
	dirichletbd = 2;
	dirichletmidbd = 3;
	
	pbcpar.we.inds = pdbcfull.w|pdbcfull.e;
	pbcpar.sn.inds = pdbcfull.s|pdbcfull.n;
	pbcpar.io.inds = filtering.p.inner.bciofull;
	
	pbcpar.we.a11.x = neumannbd;
	pbcpar.we.a11.y = dirichletbd;
	
	pbcpar.sn.a11.x = dirichletbd;
	pbcpar.sn.a11.y = neumannbd;
	
	pbcpar.io.a11.x = neumannbd;
	pbcpar.io.a11.y = dirichletbd;
	
	Lp = laplacian2(nx,ny,hx,hy,1,-1,pbcpar);
	Lp = filtering.p.inner.filterMat*Lp*filtering.p.inner.filterMat';
	Lp(1,1) = 3/2*Lp(1,1);
	%Lp(filtering.p.inner.dbc.ci,filtering.p.inner.dbc.ci) = 3/2*Lp(filtering.p.inner.dbc.ci,filtering.p.inner.dbc.ci);
	perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
	
	ubcpar.we.inds = udbcfull.w|udbcfull.e;
	ubcpar.sn.inds = udbcfull.s|udbcfull.n;
	ubcpar.io.inds = filtering.u.inner.bciofull&~udbcfull.e;
	
	ubcpar.we.a11.x = dirichletmidbd;
	ubcpar.we.a11.y = dirichletbd;
	
	ubcpar.sn.a11.x = dirichletbd;
	ubcpar.sn.a11.y = neumannbd;
	
	ubcpar.io.a11.x = dirichletmidbd;
	ubcpar.io.a11.y = dirichletbd;
	
	Lu = laplacian2(nx-1,ny,hx,hy,1,-1,ubcpar);
	Lu = par.dt/par.Re*Lu + speye(size(Lu,1));
	Lu = filtering.u.inner.filterMat*Lu*filtering.u.inner.filterMat';
	peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
	
	vbcpar.we.inds = vdbcfull.w|vdbcfull.e;
	vbcpar.sn.inds = vdbcfull.s|vdbcfull.n;
	vbcpar.io.inds = filtering.v.inner.bciofull;
	
	vbcpar.we.a11.x = neumannbd;
	vbcpar.we.a11.y = dirichletbd;
	
	vbcpar.sn.a11.x = dirichletbd;
	vbcpar.sn.a11.y = dirichletmidbd;
	
	vbcpar.io.a11.x = neumannbd;
	vbcpar.io.a11.y = dirichletbd;	
	
	Lv = laplacian2(nx,ny-1,hx,hy,1,-1,vbcpar);
	Lv = par.dt/par.Re*Lv + speye(size(Lv,1));
	Lv = filtering.v.inner.filterMat*Lv*filtering.v.inner.filterMat';
	perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
	
	qbcpar.we.inds = qdbcfull.w|qdbcfull.e;
	qbcpar.sn.inds = qdbcfull.s|qdbcfull.n;
	qbcpar.io.inds = filtering.q.inner.bciofull;
	
	qbcpar.we.a11.x = dirichletbd;
	qbcpar.we.a11.y = dirichletbd;
	
	qbcpar.sn.a11.x = dirichletbd;
	qbcpar.sn.a11.y = dirichletbd;
	
	qbcpar.io.a11.x = dirichletbd;
	qbcpar.io.a11.y = dirichletbd;
	Lq = laplacian2(nx-1,ny-1,hx,hy,1,-1,qbcpar);
	Lq = filtering.q.inner.filterMat*Lq*filtering.q.inner.filterMat';
	perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	if(par.useGPU)
		
		Rp = gpuArray(Lp);
		Rpt = Rp';
		
		Ru = gpuArray(Lu);
		Rut = Ru';
		
		Rv = gpuArray(Lv);
		Rvt = Rv';
		
		Rq = gpuArray(Lq);
		Rqt = Rq';
		
		res.Ugpu = gpuArray(grids.q.inner.Xmesh*0);
		res.Vgpu = res.Ugpu;
	end
	
	uValMatInE = [zeros([ny,1]),filtering.u.inner.valmat,zeros([ny,1])];
	uValMatInE = logical([zeros([1,nx+1]);uValMatInE;zeros([1,nx+1])]);
	
	vValMatInE = [zeros([ny-1,1]),filtering.v.inner.valmat,zeros([ny-1,1])];
	vValMatInE = logical([zeros([1,nx+2]);vValMatInE;zeros([1,nx+2])]);
	
	pValMatInE = [zeros([ny,1]),filtering.p.inner.valmat,zeros([ny,1])];
	pValMatInE = logical([zeros([1,nx+2]);pValMatInE;zeros([1,nx+2])]);
	
	qValMatInE = [zeros([ny-1,1]),filtering.q.inner.valmat,zeros([ny-1,1])];
	qValMatInE = logical([zeros([1,nx+1]);qValMatInE;zeros([1,nx+1])]);
	
% 	pOnMat = reshape(filtering.p.inner.onfull,[nx,ny])';
% 	Pselin.W = circshift(Pdbcfull.c,1,2)&pOnMat&~(grids.p.inner.Xmesh == min(grids.p.inner.xinit));
% 	Pselin.E = circshift(Pdbcfull.c,-1,2)&pOnMat&~(grids.p.inner.Xmesh == max(grids.p.inner.xinit));
% 	Pselin.S = circshift(Pdbcfull.c,1)&pOnMat&~(grids.p.inner.Ymesh == min(grids.p.inner.yinit));
% 	Pselin.N = circshift(Pdbcfull.c,-1)&pOnMat&~(grids.p.inner.Ymesh == max(grids.p.inner.yinit));
% 	Pselin.NW.c = circshift(Pselin.W,-1,2)&circshift(Pselin.N,1);
% 	Pselin.NE.c = circshift(Pselin.E,1,2)&circshift(Pselin.N,1);
% 	Pselin.SW.c = circshift(Pselin.W,-1,2)&circshift(Pselin.S,-1);
% 	Pselin.SE.c = circshift(Pselin.E,1,2)&circshift(Pselin.S,-1);
% 	Pselin.NW.N = Pselin.N&circshift(Pselin.NW.c,-1);
% 	Pselin.NW.W = Pselin.W&circshift(Pselin.NW.c,1,2);
% 	Pselin.NE.N = Pselin.N&circshift(Pselin.NE.c,-1);
% 	Pselin.NE.E = Pselin.E&circshift(Pselin.NE.c,-1,2);
% 	Pselin.SW.S = Pselin.S&circshift(Pselin.SW.c,1);
% 	Pselin.SW.W = Pselin.W&circshift(Pselin.SW.c,1,2);
% 	Pselin.SE.S = Pselin.S&circshift(Pselin.SE.c,1);
% 	Pselin.SE.E = Pselin.E&circshift(Pselin.SE.c,-1,2);
	
	Pxsel = circshift(Pdbcfull.w,-1,2)&~circshift(reshape(filtering.p.inner.bciofull,[nx,ny])',-1,2);
	Pxsel = Pxsel(:,1:end-1);
	Pysel = Pdbcfull.n|circshift(Pdbcfull.s,-1,1);
	Pysel = Pysel(1:end-1,:);
	
	Uselin.N = circshift(UdbcfullE.n,-1)&~(grids.u.outer.Ymesh == min(grids.u.outer.yinit));
	Uselin.S = circshift(UdbcfullE.s,1)&~(grids.u.outer.Ymesh == max(grids.u.outer.yinit));
	Vselin.W = circshift(VdbcfullE.w,1,2)&~(grids.v.outer.Xmesh == max(grids.v.outer.xinit));
	Vselin.E = circshift(VdbcfullE.e,-1,2)&~(grids.v.outer.Xmesh == min(grids.v.outer.xinit));
	
	uvxsel = reshape(filtering.u.outer.valind,[nx+1,ny+2])';
	uvysel = reshape(filtering.v.outer.valind,[nx+2,ny+1])';
	
	u = zeros(size(Lu,1),1);
	v = zeros(size(Lv,1),1);
	p = zeros(size(Lp,1),1);
	q = zeros(size(Lq,1),1);
	
	numfps = par.numfps;
	pfps = zeros(numfps,1);
	cfps = zeros(numfps,1);
	cfpsi = 1;
	pfpsi = 1;
	
	pt = tic;
	ct = tic;
	
	fprintf('\n')
	lents = 0;
	for j=1:par.timesteps
		
		% nonlinear terms/explicit convection
		gamma = min(1.2*par.dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
		
		Ue(uValMatInE) = U(filtering.u.inner.valmat);
		Ue(uselEw|uselEe) = UbcE(uselEw|uselEe);
		Ue(UdbcfullE.n|UdbcfullE.s) = 2*UbcE(UdbcfullE.n|UdbcfullE.s)-Ue(Uselin.N|Uselin.S);
		
		Ve(vValMatInE) = V(filtering.v.inner.valmat);
		Ve(vselEs|vselEn) = VbcE(vselEs|vselEn);
		Ve(VdbcfullE.w|VdbcfullE.e) = 2*VbcE(VdbcfullE.w|VdbcfullE.e)-Ve(Vselin.W|Vselin.E);
		
		Ua = mvgavg(Ue);
		Ud = diff(Ue)/2;
		Va = mvgavg(Ve,2);
		Vd = diff(Ve')'/2;
		UVx = diff((Ua.*Va-gamma*abs(Ua).*Vd)')'/hx;
		UVy = diff((Ua.*Va-gamma*Ud.*abs(Va)))/hy;
		
		UVx(VdbcfullE.w(:,2:end-1)|VdbcfullE.e(:,2:end-1)|~uvysel(:,2:end-1)) = 0;
		UVy(UdbcfullE.s(2:end-1,:)|UdbcfullE.n(2:end-1,:)|~uvxsel(2:end-1,:)) = 0;
		
		Ue2 = Ue;
		Ue2(~uValMatInE&(UdbcfullE.n|UdbcfullE.s)) = 0;
		Ue2 = Ue2(2:end-1,:);
		Ve2 = Ve;
		Ve2(~vValMatInE&(VdbcfullE.w|VdbcfullE.e)) = 0;
		Ve2 = Ve2(:,2:end-1);
		Ua = mvgavg(Ue2,2);
		Ud = diff(Ue2')'/2;
		Va = mvgavg(Ve2);
		Vd = diff(Ve2)/2;
		U2x = diff((Ua.^2-gamma*abs(Ua).*Ud)')'/hx;
		V2y = diff((Va.^2-gamma*abs(Va).*Vd))/hy;
		
		U = U-par.dt*(UVy(:,2:end-1)+U2x);
		V = V-par.dt*(UVx(2:end-1,:)+V2y);
		
		%U(Udbcfull.io) = Ubc(Udbcfull.io);
		%V(Vdbcfull.io) = Vbc(Vdbcfull.io);
		
		% implicit diffusion
		rhsu = reshape((U+Ubc)',[],1);
		rhsu = rhsu(filtering.u.inner.valind);
		u(peru) = (Ru\(Rut\rhsu(peru)));
		U = reshape(filtering.u.inner.filterMat'*u,nx-1,ny)';
		
		rhsv = reshape((V+Vbc)',[],1);
		rhsv = rhsv(filtering.v.inner.valind);
		v(perv) = (Rv\(Rvt\rhsv(perv)));
		V = reshape(filtering.v.inner.filterMat'*v,nx,ny-1)';
		
		% pressure correction
		Uep = Ue;
		Uep(uValMatInE) = U(filtering.u.inner.valmat);
		Ue(uselEw|uselEe) = UbcE(uselEw|uselEe);
		Uep(UdbcfullE.s|UdbcfullE.n) = 0;
		Uep = Uep(2:end-1,:);
		
		Vep = Ve;
		Vep(vValMatInE) = V(filtering.v.inner.valmat);
		Vep(vselEs|vselEn) = VbcE(vselEs|vselEn);
		Vep(VdbcfullE.w|VdbcfullE.e) = 0;
		Vep = Vep(:,2:end-1);
		
		rhsp = reshape((diff(Uep')'/hx+diff(Vep)/hy)',[],1);
		rhsp = 1/par.dt*rhsp(filtering.p.inner.valind);
		p(perp) = -(Rp\(Rpt\rhsp(perp)));
		P = reshape(filtering.p.inner.filterMat'*p,nx,ny)';
% 		P(Pselin.NW.c) = (P(Pselin.NW.N)+P(Pselin.NW.W))/2;
% 		P(Pselin.NE.c) = (P(Pselin.NE.N)+P(Pselin.NE.E))/2;
% 		P(Pselin.SW.c) = (P(Pselin.SW.S)+P(Pselin.SW.W))/2;
% 		P(Pselin.SE.c) = (P(Pselin.SE.S)+P(Pselin.SE.E))/2;
% 		P(Pselin.NW.c) = 2*P(Pselin.NW.N) - P(circshift(Pselin.NW.N,-1));
% 		P(Pselin.SW.c) = 2*P(Pselin.SW.S) - P(circshift(Pselin.SW.S,1));
		
		Px = diff(P')'/hx;
		Px(Pxsel) = 0;
		Py = diff(P)/hy;
		Py(Pysel) = 0;
		U = U-par.dt*Px;
		V = V-par.dt*Py;
		
		Ue(uValMatInE) = U(filtering.u.inner.valmat);
		Ve(vValMatInE) = V(filtering.v.inner.valmat);
		
		for i=1:numel(par.wesn)
			Ue(UdbcfullE.(par.wesn{i})) = UbcE(UdbcfullE.(par.wesn{i}));
			Ve(VdbcfullE.(par.wesn{i})) = VbcE(VdbcfullE.(par.wesn{i}));
		end
		
		if((j==1 || mod(j,par.plotoniter)==0)&&~par.noplot)
			
			%stream function
% 			rhsq = reshape((diff(U)/hy-diff(V')'/hx)',[],1);
% 			%rhsq(filtering.q.inner.onfull) = rhs.q.inner(filtering.q.inner.on);
% 			rhsq = rhsq(filtering.q.inner.valind);
% 			q(perq) = Rq\(Rqt\rhsq(perq));
% 			Q = reshape(filtering.q.inner.filterMat'*q,nx-1,ny-1)';
			
			[res.Ue,res.Ve] = makeQuiverData(Ue,Ve,QdbcfullE,QUbcE,QVbcE,reshape(filtering.q.outer.valind,[nx+1,ny+1])');
			
			if(par.useGPU)
				res.Ugpu(filtering.q.inner.valmat) = res.Ue(qValMatInE);
				res.Vgpu(filtering.q.inner.valmat) = res.Ve(qValMatInE);
				
				%a = res.Ue;
				%b = res.U;
				res.P = gather(P);
				Qu = cumsum(res.Ugpu*hy);
				%Qv = -cumsum(res.V'*hx)';
				res.Q = Qu;
				
				res.Qe = cumsum(res.Ue*hy);
				
				res.U = gather(res.Ugpu);
				res.V = gather(res.Vgpu);
				res.P = gather(res.P);
				res.Q = gather(res.Q);
				res.Ue = gather(res.Ue);
				res.Ve = gather(res.Ve);
				res.Qe = gather(res.Qe);
			else
				
				res.U = zeros(ny-1,nx-1);
				res.V = res.U;
				res.U(filtering.q.inner.valmat) = res.Ue(qValMatInE);
				res.V(filtering.q.inner.valmat) = res.Ve(qValMatInE);
				
				%a = res.Ue;
				%b = res.U;
				res.P = P;
				Qu = cumsum(res.U*hy);
				Qv = -cumsum(res.V'*hx)';
				res.Q = Qu;
				%res.Q = Q;
				
				res.Qe = cumsum(res.Ue*hy);
			end
			
			if(par.useinterp)

				res.U  = interp2(grids.q.inner.X,grids.q.inner.Y,res.U ,grids.plot.q.inner.X,grids.plot.q.inner.Y);
				res.V  = interp2(grids.q.inner.X,grids.q.inner.Y,res.V ,grids.plot.q.inner.X,grids.plot.q.inner.Y);
				res.P  = interp2(grids.p.inner.X,grids.p.inner.Y,res.P ,grids.plot.p.inner.X,grids.plot.p.inner.Y);
				res.Q  = interp2(grids.q.inner.X,grids.q.inner.Y,res.Q ,grids.plot.q.inner.X,grids.plot.q.inner.Y);
				res.Ue = interp2(grids.q.outer.X,grids.q.outer.Y,res.Ue,grids.plot.q.outer.X,grids.plot.q.outer.Y);
				res.Ve = interp2(grids.q.outer.X,grids.q.outer.Y,res.Ve,grids.plot.q.outer.X,grids.plot.q.outer.Y);
				res.Qe = interp2(grids.q.outer.X,grids.q.outer.Y,res.Qe,grids.plot.q.outer.X,grids.plot.q.outer.Y);

				PlotNS(grids.plot,filtering.plot,res,par)

			else

				PlotNS(grids,filtering,res,par);

			end
			
			pfps(pfpsi) = 1/toc(pt);
			pt = tic;
			if(pfpsi == numfps)
				pfpsi = 1;
			else
				pfpsi = pfpsi+1;
			end
			
		end
		
		cfps(cfpsi) = 1/toc(ct);
		ct = tic;
		if(cfpsi == numfps)
			cfpsi = 1;
		else
			cfpsi = cfpsi+1;
		end	
		
		timestr = sprintf('Time step: %5.2f/%5.2fs | Plotted FPS: %d fps | Computed FPS: %d fps',j*par.dt,par.tf,round(sum(pfps)/(numfps+min(j-numfps,0))),round(sum(cfps)/(numfps+min(j-numfps,0))));
		fprintf([repmat('\b',1,lents) timestr]);
		lents = numel(timestr);
		
	end
	fprintf('\n');
end

function M = gridify(v,nx,ny)
	M = reshape(v,[nx,ny])';
end

function v = vectify(M)
	v = reshape(M',[],1);
end

function [Unew,Vnew] = makeQuiverData(U,V,QdbcfullE,QUbcE,QVbcE,qvalmat)
	Unew = mvgavg(U);
	Vnew = mvgavg(V,2);
	Unew(QdbcfullE.s|QdbcfullE.n|QdbcfullE.c) = QUbcE(QdbcfullE.s|QdbcfullE.n|QdbcfullE.c);
	Vnew(QdbcfullE.w|QdbcfullE.e|QdbcfullE.c) = QVbcE(QdbcfullE.w|QdbcfullE.e|QdbcfullE.c);
	Unew(~(qvalmat)) = 0;
	Vnew(~(qvalmat)) = 0;
	%Len = sqrt(Unew.^2+Vnew.^2+eps);
	%Unew = Unew./Len;
	%Vnew = Vnew./Len;
end

