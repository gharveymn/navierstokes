function [figs,mat,vec] = NSPrim(par,grids,filtering,rhs,figs)
	%NSPRIM primitive formulation of NS
	
	nx = grids.p.nx;
	ny = grids.p.ny;
	h = par.h;
	
	U = 0*grids.u.inner.Xmesh;
	V = 0*grids.v.inner.Ymesh;
	Ue = 0*grids.u.outer.Xmesh;
	Ve = 0*grids.v.outer.Ymesh;
	
	Ubc = par.dt/(par.h*par.Re)*reshape(filtering.u.inner.filterMat'*rhs.inner.u,[nx-1,ny])';
	Vbc = par.dt/(par.h*par.Re)*reshape(filtering.v.inner.filterMat'*rhs.inner.v*0,[nx,ny-1])';
	
	Ubcfull = reshape(filtering.u.outer.filterMat'*rhs.outer.u,[nx+1,ny+2])';
	Vbcfull = reshape(filtering.v.outer.filterMat'*rhs.outer.v,[nx+2,ny+1])';
	
	pdbcfull = filtering.p.inner.dbcfull;
	Lp = laplacian2(nx,ny,h,1,1,1,1,pdbcfull.w|pdbcfull.e,pdbcfull.s|pdbcfull.n,filtering.p.inner.bcio,-1);
	Lp = filtering.p.inner.filterMat*Lp*filtering.p.inner.filterMat';
	Lp(1,1) = 3/2*Lp(1,1);
	%perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
	
	udbcfull = filtering.u.inner.dbcfull;
	Lu = laplacian2(nx-1,ny,h,2,3,2,udbcfull.w|udbcfull.e,udbcfull.s|udbcfull.n,filtering.u.inner.bcio,-1);
	Lu = par.dt/par.Re*Lu + speye(size(Lu,1));
	%peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
	
	vdbcfull = filtering.v.inner.dbcfull;
	Lv = laplacian2(nx,ny-1,h,3,2,3,vdbcfull.w|vdbcfull.e,vdbcfull.s|vdbcfull.n,filtering.v.inner.bcio,-1);
	Lv = par.dt/par.Re*Lv + speye(size(Lv,1));
	%perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
	
	qdbcfull = filtering.q.inner.dbcfull;
	Lq = laplacian2(nx-1,ny-1,h,2,2,2,qdbcfull.w|qdbcfull.e,qdbcfull.s|qdbcfull.n,filtering.q.inner.bcio,-1);
	%perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	for i=1:par.timesteps
		
		%nonlinear terms
		gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
		Ue(2:end-1,2:end-1) = U;
		udbcfull = filtering.u.outer.dbcfull;
		U(udbcfull.w) = 0;
		U(udbcfull.e) = 0;
		U(udbcfull.n) = 0-U(circshift(udbcfull.n,-1)&~(udbcfull.w|udbcfull.e|udbcfull.s));
		U(udbcfull.s) = 0-U(circshift(udbcfull.s,1)&~(udbcfull.w|udbcfull.e|udbcfull.n));
		sel = reshape(filtering.u.outer.filterMat'*filtering.u.outer.bcio,[nx+1,ny+2]);
		U(sel) = Ubcfull(sel);
		
		V(vdbcfull.s) = 0;
		V(vdbcfull.n) = 0;
		V(vdbcfull.w) = 0-V(circshift(vdbcfull.n,-nx)&~(vdbcfull.s|vdbcfull.n|vdbcfull.e));
		V(vdbcfull.e) = 0-V(circshift(vdbcfull.s,nx)&~(vdbcfull.s|vdbcfull.n|vdbcfull.w));
		sel = reshape(filtering.v.outer.filterMat'*filtering.v.outer.bcio,[nx+1,ny+2]);
		V(sel) = Vbcfull(sel);	
		
		Ua = mvgavg(Ue);
		Ud = diff(Ue/2);
		Va = mvgavg(Ve,2);
		Vd = diff(Ve/2);
		UVx = diff((Ua.*Va-gamma*abs(Ua).*Vd)')'/h;
		UVy = diff((Ua.*Va-gamma*Ud.*abs(Va)))/h;
		
		Ua = avg(Ue(:,2:end-1));   Ud = diff(Ue(:,2:end-1))/2;
		Va = avg(Ve(2:end-1,:)')'; Vd = diff(Ve(2:end-1,:)')'/2;
		U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
		V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
		U = U-dt*(UVy(2:end-1,:)+U2x);
		V = V-dt*(UVx(:,2:end-1)+V2y);	
		
		if(~exist('figs','var'))
			[figs,mat,vec] = InPost(qnew,bc,grids,filtering,par);
		else
			[figs,mat,vec] = InPost(qnew,bc,grids,filtering,par,figs);
		end
	end
	
	
end

