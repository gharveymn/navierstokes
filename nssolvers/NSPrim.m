function [grids,filtering,res,par,figs] = NSPrim(par,grids,filtering,rhs,figs)
	%NSPRIM primitive formulation of NS
	
	nx = grids.p.nx;
	ny = grids.p.ny;
	h = par.h;
	
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
	%gridify
	Udbcfull.w = gridify(udbcfull.w,nx-1,ny);
	Udbcfull.e = gridify(udbcfull.e,nx-1,ny);
	Udbcfull.s = gridify(udbcfull.s,nx-1,ny);
	Udbcfull.n = gridify(udbcfull.n,nx-1,ny);
	Udbcfull.c = gridify(udbcfull.c,nx-1,ny);
	
	Vdbcfull.w = gridify(vdbcfull.w,nx,ny-1);
	Vdbcfull.e = gridify(vdbcfull.e,nx,ny-1);
	Vdbcfull.s = gridify(vdbcfull.s,nx,ny-1);
	Vdbcfull.n = gridify(vdbcfull.n,nx,ny-1);
	Vdbcfull.c = gridify(vdbcfull.c,nx,ny-1);
	
	Pdbcfull.w = gridify(pdbcfull.w,nx,ny);
	Pdbcfull.e = gridify(pdbcfull.e,nx,ny);
	Pdbcfull.s = gridify(pdbcfull.s,nx,ny);
	Pdbcfull.n = gridify(pdbcfull.n,nx,ny);
	Pdbcfull.c = gridify(pdbcfull.c,nx,ny);
	
	%outer
	udbcfullE = filtering.u.outer.dbcfull;
	vdbcfullE = filtering.v.outer.dbcfull;
	%gridify
	UdbcfullE.w = gridify(udbcfullE.w,nx+1,ny+2);
	UdbcfullE.e = gridify(udbcfullE.e,nx+1,ny+2);
	UdbcfullE.s = gridify(udbcfullE.s,nx+1,ny+2);
	UdbcfullE.n = gridify(udbcfullE.n,nx+1,ny+2);
	UdbcfullE.c = gridify(udbcfullE.c,nx+1,ny+2);
	
	VdbcfullE.w = gridify(vdbcfullE.w,nx+2,ny+1);
	VdbcfullE.e = gridify(vdbcfullE.e,nx+2,ny+1);
	VdbcfullE.s = gridify(vdbcfullE.s,nx+2,ny+1);
	VdbcfullE.n = gridify(vdbcfullE.n,nx+2,ny+1);
	VdbcfullE.c = gridify(vdbcfullE.c,nx+2,ny+1);
	
	uselEw = UdbcfullE.w&~(UdbcfullE.s|UdbcfullE.n);
	uselEe = UdbcfullE.e&~(UdbcfullE.s|UdbcfullE.n);
	%uselEs = UdbcfullE.s&~(UdbcfullE.w|UdbcfullE.e);
	%uselEn = UdbcfullE.n&~(UdbcfullE.w|UdbcfullE.e);	
	
	%vselEw = VdbcfullE.w&~(VdbcfullE.s|VdbcfullE.n);
	%vselEe = VdbcfullE.e&~(VdbcfullE.s|VdbcfullE.n);
	vselEs = VdbcfullE.s&~(VdbcfullE.w|VdbcfullE.e);
	vselEn = VdbcfullE.n&~(VdbcfullE.w|VdbcfullE.e);
	
	%inner boundary cond
	Urhs = reshape(filtering.u.inner.filterMat'*rhs.inner.u,[nx-1,ny])';
	Ubc = Urhs;
	Ubc(Udbcfull.n|Udbcfull.s) = 2*par.dt/(par.h^2*par.Re)*Urhs(Udbcfull.n|Udbcfull.s);
	Ubc(Udbcfull.w|Udbcfull.e) = par.dt/(par.h^2*par.Re)*Urhs(Udbcfull.w|Udbcfull.e);
	
	Vrhs = reshape(filtering.v.inner.filterMat'*rhs.inner.v*0,[nx,ny-1])';
	Vbc = Vrhs;
	Vbc(Vdbcfull.w|Vdbcfull.e) = 2*par.dt/(par.h^2*par.Re)*Vrhs(Vdbcfull.w|Vdbcfull.e);
	Vbc(Vdbcfull.n|Vdbcfull.s) = par.dt/(par.h^2*par.Re)*Vrhs(Vdbcfull.s|Vdbcfull.n);
	
	%outer boundary cond
	Ubcfull = reshape(filtering.u.outer.filterMat'*rhs.outer.u,[nx+1,ny+2])';
	Vbcfull = reshape(filtering.v.outer.filterMat'*rhs.outer.v*0,[nx+2,ny+1])';
	
	Lp = laplacian2(nx,ny,h,1,1,1,1,pdbcfull.w|pdbcfull.e,pdbcfull.s|pdbcfull.n,filtering.p.inner.bciofull,-1);
	Lp = filtering.p.inner.filterMat*Lp*filtering.p.inner.filterMat';
	Lp(1,1) = 3/2*Lp(1,1);
	%perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
	
	Lu = laplacian2(nx-1,ny,h,3,2,3,1,udbcfull.w|udbcfull.e,udbcfull.s|udbcfull.n,filtering.u.inner.bciofull,-1);
	Lu = par.dt/par.Re*Lu + speye(size(Lu,1));
	Lu = filtering.u.inner.filterMat*Lu*filtering.u.inner.filterMat';
	%peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
	
	Lv = laplacian2(nx,ny-1,h,2,3,2,1,vdbcfull.w|vdbcfull.e,vdbcfull.s|vdbcfull.n,filtering.v.inner.bciofull,-1);
	Lv = par.dt/par.Re*Lv + speye(size(Lv,1));
	Lv = filtering.v.inner.filterMat*Lv*filtering.v.inner.filterMat';
	%perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
	
	Lq = laplacian2(nx-1,ny-1,h,2,2,2,1,qdbcfull.w|qdbcfull.e,qdbcfull.s|qdbcfull.n,filtering.q.inner.bciofull,-1);
	Lq = filtering.q.inner.filterMat*Lq*filtering.q.inner.filterMat';
	%perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	if(par.useGPU)
		Lp = gpuArray(Lp);
		Lu = gpuArray(Lu);
		Lv = gpuArray(Lv);
		Lq = gpuArray(Lq);
	end
	
	uValMat = reshape(filtering.u.inner.valind,[nx-1,ny])';
	uValMatE = [zeros([ny,1]),uValMat,zeros([ny,1])];
	uValMatE = logical([zeros([1,nx+1]);uValMatE;zeros([1,nx+1])]);
	uvalindE = reshape(uValMatE',[],1);
	
	vValMat = reshape(filtering.v.inner.valind,[nx,ny-1])';
	vValMatE = [zeros([ny-1,1]),vValMat,zeros([ny-1,1])];
	vValMatE = logical([zeros([1,nx+2]);vValMatE;zeros([1,nx+2])]);
	vvalindE = reshape(vValMatE',[],1);
	
	pValMat = reshape(filtering.p.inner.valind,[nx,ny])';
	pValMatE = [zeros([ny,1]),pValMat,zeros([ny,1])];
	pValMatE = logical([zeros([1,nx+2]);pValMatE;zeros([1,nx+2])]);
	pvalindE = reshape(pValMatE',[],1);
	
	qValMat = reshape(filtering.q.inner.valind,[nx-1,ny-1])';
	qValMatE = [zeros([ny-1,1]),qValMat,zeros([ny-1,1])];
	qValMatE = logical([zeros([1,nx+1]);qValMatE;zeros([1,nx+1])]);
	qvalindE = reshape(qValMatE',[],1);
	
	Pxsel = circshift(Pdbcfull.w,-1,2)&~circshift(reshape(filtering.p.inner.bciofull,[nx,ny])',-1,2);
	Pxsel = Pxsel(:,1:end-1);
	Pysel = Pdbcfull.n|circshift(Pdbcfull.s,-1,1);
	Pysel = Pysel(1:end-1,:);
	
	UselinN = circshift(UdbcfullE.n,-1)&~(grids.u.outer.Ymesh == min(grids.u.outer.yinit));
	UselinS = circshift(UdbcfullE.s,1)&~(grids.u.outer.Ymesh == max(grids.u.outer.yinit));
	VselinW = circshift(VdbcfullE.w,1,2)&~(grids.v.outer.Xmesh == max(grids.v.outer.xinit));
	VselinE = circshift(VdbcfullE.e,-1,2)&~(grids.v.outer.Xmesh == min(grids.v.outer.xinit));
	
	fprintf('\n')
	lents = 0;
	for i=1:par.timesteps
		
		%nonlinear terms
		gamma = min(1.2*par.dt*max(max(max(abs(U)))/par.h,max(max(abs(V)))/par.h),1);
		
		Ue(uValMatE) = U(uValMat);
		Ue(uselEw) = Ubcfull(uselEw);
		Ue(uselEe) = Ubcfull(uselEe);
		Ue(UdbcfullE.n) = 2*Ubcfull(UdbcfullE.n)-Ue(UselinN);
		Ue(UdbcfullE.s) = 2*Ubcfull(UdbcfullE.s)-Ue(UselinS);
		
		Ve(vValMatE) = V(vValMat);
		Ve(vselEs) = Vbcfull(vselEs);
		Ve(vselEn) = Vbcfull(vselEn);
		Ve(VdbcfullE.w) = 2*Vbcfull(VdbcfullE.w)-Ve(VselinW);
		Ve(VdbcfullE.e) = 2*Vbcfull(VdbcfullE.e)-Ve(VselinE);
		
		Ua = mvgavg(Ue);
		Ud = diff(Ue)/2;
		Va = mvgavg(Ve,2);
		Vd = diff(Ve')'/2;
		UVx = diff((Ua.*Va-gamma*abs(Ua).*Vd)')'/h;
		UVy = diff((Ua.*Va-gamma*Ud.*abs(Va)))/h;
		
		Ue2 = Ue;
		Ue2(~uValMatE&(UdbcfullE.n|UdbcfullE.s)) = 0;
		Ue2 = Ue2(2:end-1,:);
		Ve2 = Ve;
		Ve2(~vValMatE&(VdbcfullE.w|VdbcfullE.e)) = 0;
		Ve2 = Ve2(:,2:end-1);
		Ua = mvgavg(Ue2,2);
		Ud = diff(Ue2')'/2;
		Va = mvgavg(Ve2);
		Vd = diff(Ve2)/2;
		U2x = diff((Ua.^2-gamma*abs(Ua).*Ud)')'/par.h;
		V2y = diff((Va.^2-gamma*abs(Va).*Vd))/par.h;
		
		U = U-par.dt*(UVy(:,2:end-1)+U2x);
		V = V-par.dt*(UVx(2:end-1,:)+V2y);
		
		% implicit viscosity
		rhs = reshape((U+Ubc)',[],1);
		rhs = rhs(filtering.u.inner.valind);
		u = Lu\rhs;
		U = reshape((1*filtering.u.inner.filterMat')*u,nx-1,ny)';
		
		rhs = reshape((V+Vbc)',[],1);
		rhs = rhs(filtering.v.inner.valind);
		v = Lv\rhs;
		V = reshape((1*filtering.v.inner.filterMat')*v,nx,ny-1)';
		
		%pressure correction
		Uep = Ue;
		Uep(uValMatE) = U(uValMat);
		Uep(UdbcfullE.s|UdbcfullE.n) = 0;
		Uep = Uep(2:end-1,:);
		
		Vep = Ve;
		Vep(vValMatE) = V(vValMat);					
		Vep(VdbcfullE.w|VdbcfullE.e) = 0;
		Vep = Vep(:,2:end-1);
		
		rhs = reshape((diff(Uep')'+diff(Vep))'/par.h,[],1);
		p = -Lp\rhs(filtering.p.inner.valind);
		P = reshape((1*filtering.p.inner.filterMat')*p,nx,ny)';
		Px = diff(P')'/par.h;
		%Px(Pxsel) = 0;
		Py = diff(P)/par.h;
		%Py(Pysel) = 0;
		U = U-Px;
		V = V-Py;
		
		%stream function
		rhs = reshape((diff(U)-diff(V')')'/par.h,[],1);
		q = Lq\rhs(filtering.q.inner.valind);
		Q = reshape((1*filtering.q.inner.filterMat')*q,nx-1,ny-1)';
		
		res.U = gather(U);
		res.V = gather(V);
		res.P = gather(P);
		res.Q = gather(Q);
		
		if(i==1 || mod(i,par.plotoniter)==0)
			if(~exist('figs','var'))
				figs = PlotNS(grids,filtering,res,par);
			else
				figs = PlotNS(grids,filtering,res,par,figs);
			end
		end
		
		
		
		timestr = ['Time step: ' num2str(i*par.dt,'%5.2f') '/' num2str(par.tf,'%5.2f') 's'];
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

