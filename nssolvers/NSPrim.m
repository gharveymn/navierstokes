function [figs,mat,vec] = NSPrim(par,grids,filtering,rhs,figs)
	%NSPRIM primitive formulation of NS
	
	nx = grids.p.nx;
	ny = grids.p.ny;
	h = par.h;
	
	U = zeros([ny,nx-1]);
	V = zeros([ny-1,nx]);
	Ue = zeros([ny+2,nx+1]);
	Ve = zeros([ny+1,nx+2]);
	
	Ubc = par.dt/(par.h*par.Re)*reshape(filtering.u.inner.filterMat'*rhs.inner.u,[nx-1,ny])';
	Vbc = par.dt/(par.h*par.Re)*reshape(filtering.v.inner.filterMat'*rhs.inner.v*0,[nx,ny-1])';
	
	Ubcfull = reshape(filtering.u.outer.filterMat'*rhs.outer.u,[nx+1,ny+2])';
	Vbcfull = reshape(filtering.v.outer.filterMat'*rhs.outer.v*0,[nx+2,ny+1])';
	
	pdbcfull = filtering.p.inner.dbcfull;
	Lp = laplacian2(nx,ny,h,1,1,1,1,pdbcfull.w|pdbcfull.e,pdbcfull.s|pdbcfull.n,filtering.p.inner.bciofull,-1);
	Lp = filtering.p.inner.filterMat*Lp*filtering.p.inner.filterMat';
	Lp(1,1) = 3/2*Lp(1,1);
	%perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
	
	udbcfull = filtering.u.inner.dbcfull;
	Lu = laplacian2(nx-1,ny,h,3,2,3,1,udbcfull.w|udbcfull.e,udbcfull.s|udbcfull.n,filtering.u.inner.bciofull,-1);
	Lu = par.dt/par.Re*Lu + speye(size(Lu,1));
	Lu = filtering.u.inner.filterMat*Lu*filtering.u.inner.filterMat';
	%peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
	
	vdbcfull = filtering.v.inner.dbcfull;
	Lv = laplacian2(nx,ny-1,h,2,3,2,1,vdbcfull.w|vdbcfull.e,vdbcfull.s|vdbcfull.n,filtering.v.inner.bciofull,-1);
	Lv = par.dt/par.Re*Lv + speye(size(Lv,1));
	Lv = filtering.v.inner.filterMat*Lv*filtering.v.inner.filterMat';
	%perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
	
	qdbcfull = filtering.q.inner.dbcfull;
	Lq = laplacian2(nx-1,ny-1,h,2,2,2,1,qdbcfull.w|qdbcfull.e,qdbcfull.s|qdbcfull.n,filtering.q.inner.bciofull,-1);
	Lq = filtering.q.inner.filterMat*Lq*filtering.q.inner.filterMat';
	%perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
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
	
	udbcfullE = filtering.u.outer.dbcfull;
	vdbcfullE = filtering.v.outer.dbcfull;
	%gridify
	udbcfullE.w = gridify(udbcfullE.w,nx+1,ny+2);
	udbcfullE.e = gridify(udbcfullE.e,nx+1,ny+2);
	udbcfullE.s = gridify(udbcfullE.s,nx+1,ny+2);
	udbcfullE.n = gridify(udbcfullE.n,nx+1,ny+2);
	udbcfullE.c = gridify(udbcfullE.c,nx+1,ny+2);
	
	vdbcfullE.w = gridify(vdbcfullE.w,nx+2,ny+1);
	vdbcfullE.e = gridify(vdbcfullE.e,nx+2,ny+1);
	vdbcfullE.s = gridify(vdbcfullE.s,nx+2,ny+1);
	vdbcfullE.n = gridify(vdbcfullE.n,nx+2,ny+1);
	vdbcfullE.c = gridify(vdbcfullE.c,nx+2,ny+1);
	
	for i=1:par.timesteps
		
		%nonlinear terms
		gamma = min(1.2*par.dt*max(max(max(abs(U)))/par.h,max(max(abs(V)))/par.h),1);
		
		Ue(uValMatE) = U(uValMat);
		Ue(udbcfullE.w) = Ubcfull(udbcfullE.w);
		Ue(udbcfullE.e) = Ubcfull(udbcfullE.e);
		Ue(udbcfullE.n) = 2*Ubcfull(udbcfullE.n)-Ue(udbcfullE.n);
		Ue(udbcfullE.s) = 2*Ubcfull(udbcfullE.s)-Ue(udbcfullE.s);
		sel = logical(reshape(filtering.u.outer.filterMat'*(1*filtering.u.outer.bcio),[nx+1,ny+2])');
		Ue(sel) = Ubcfull(sel);
		
		Ve(vValMatE) = V(vValMat);
		Ve(vdbcfullE.s) = Vbcfull(vdbcfullE.s);
		Ve(vdbcfullE.n) = Vbcfull(vdbcfullE.n);
		Ve(vdbcfullE.w) = 2*Vbcfull(vdbcfullE.w)-Ve(vdbcfullE.w);
		Ve(vdbcfullE.e) = 2*Vbcfull(vdbcfullE.e)-Ve(vdbcfullE.e);
		sel = logical(reshape(filtering.v.outer.filterMat'*(1*filtering.v.outer.bcio),[nx+2,ny+1])');
		Ve(sel) = Vbcfull(sel);
		
		Ua = mvgavg(Ue);
		Ud = diff(Ue)/2;
		Va = mvgavg(Ve,2);
		Vd = diff(Ve')'/2;
		UVx = diff((Ua.*Va-gamma*abs(Ua).*Vd)')'/h;
		UVy = diff((Ua.*Va-gamma*Ud.*abs(Va)))/h;
		
		Ue2 = Ue;
		Ue2(~uValMatE&(udbcfullE.n|udbcfullE.s)) = 0;
		Ue2 = Ue2(2:end-1,:);
		Ve2 = Ve;
		Ve2(~vValMatE&(vdbcfullE.w|vdbcfullE.e)) = 0;
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
		u = Lu\rhs(filtering.u.inner.valind);
		U = reshape(filtering.u.inner.filterMat'*u,nx-1,ny)';
		
		rhs = reshape((V+Vbc)',[],1);
		v = Lv\rhs(filtering.v.inner.valind);
		V = reshape(filtering.v.inner.filterMat'*v,nx,ny-1)';
		
		%pressure correction
		rhs = reshape((diff(Ue2')'+diff(Ve2))'/par.h,[],1);
		p = -Lp\rhs(filtering.p.inner.valind);
		P = reshape(filtering.p.inner.filterMat'*p,nx,ny)';
		U = U-diff(P')'/par.h;
		V = V-diff(P)/par.h;
		
		%stream function
		rhs = reshape((diff(U)-diff(V')')'/par.h,[],1);
		q = Lq\rhs(filtering.q.inner.valind);
		Q = reshape(filtering.q.inner.filterMat'*q,nx-1,ny-1)';
		
		res.U = U;
		res.V = V;
		res.P = P;
		res.Q = Q;
		
% 		figure(1)
% 		surf(grids.u.inner.Xmesh,grids.u.inner.Ymesh,U);
% 		figure(2)
% 		surf(grids.v.inner.Xmesh,grids.v.inner.Ymesh,V);
% 		figure(3)
% 		surf(grids.p.inner.Xmesh,grids.p.inner.Ymesh,P);
% 		figure(4)
% 		surf(grids.q.inner.Xmesh,grids.q.inner.Ymesh,Q);
% 		figure(5)
% 		U3 = mvgavg(U);
% 		V3 = mvgavg(V,2);
% 		Len = sqrt(U3.^2+V3.^2+eps);
% 		quiver(grids.q.inner.Xmesh,grids.q.inner.Ymesh,U3./Len,V3./Len);
% 		
% 		figure(6)
% 		contour(grids.q.inner.Xmesh,grids.q.inner.Ymesh,Q);
% 		drawnow;
		
 		if(~exist('figs','var'))
 			figs = PlotNS(grids,filtering,res,par);
 		else
 			figs = PlotNS(grids,filtering,res,par,figs);
 		end
	end
	
	
end

function M = gridify(v,nx,ny)
	M = reshape(v,[nx,ny])';
end

function v = vectify(M)
	v = reshape(M',[],1);
end

