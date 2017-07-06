function [figs,mat,vec] = NSPrim(par,grids,filtering,rhs,figs)
	%NSPRIM primitive formulation of NS
	
	nx = grids.p.nx;
	ny = grids.p.ny;
	h = par.h;
	
	U = zeros(ny,nx-1);
	V = zeros(ny-1,nx);
	
	
	
	
	Lp = laplacian2(nx,ny,h,1,1,1,filtering.p.bcw|filtering.p.bce,filtering.p.bcs|filtering.p.bcn,filtering.p.bcio,-1);
	Lp = filtering.p.filterMat*Lp*filtering.p.filterMat';
	Lp(1,1) = 3/2*Lp(1,1);
	perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
	
	Lu = laplacian2(nx-1,ny,h,2,3,2,filtering.u.bcw|filtering.u.bce,filtering.u.bcs|filtering.u.bcn,filtering.u.bcio,-1);
	Lu = par.dt/par.Re*Lu + speye(size(Lu,1));
	peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
	
	Lv = laplacian2(nx,ny-1,h,3,2,3,filtering.v.bcw|filtering.v.bce,filtering.v.bcs|filtering.v.bcn,filtering.v.bcio,-1);
	Lv = par.dt/par.Re*Lv + speye(size(Lv,1));
	perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
	
	Lq = laplacian2(nx-1,ny-1,h,2,2,2,filtering.q.bcw|filtering.q.bce,filtering.q.bcs|filtering.q.bcn,filtering.q.bcio,-1);
	perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	for i=1:par.timesteps
		
		
		
		
		if(~exist('figs','var'))
			[figs,mat,vec] = InPost(qnew,bc,grids,filtering,par);
		else
			[figs,mat,vec] = InPost(qnew,bc,grids,filtering,par,figs);
		end
	end
	
	
end

