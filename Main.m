if(~exist('testrun','var') || ~testrun)
	par = Parameters;
end


if(par.usestagger)
	
	[grids,filtering,par] = par.gridmaker(par);
	
	rhs = [];
	
	%inner
	rhs.inner.p = par.rhfunc(grids.p.inner.xmesh,grids.p.inner.ymesh);
	[rhs.inner.p,filtering.p.inner.bcio] = par.bcfunc(grids.p.inner,filtering.p.inner,rhs.inner.p,par,'p');
	filtering.p.inner.bciofull = logical(filtering.p.inner.filterMat'*(1*filtering.p.inner.bcio));
	
	rhs.inner.u = par.rhfunc(grids.u.inner.xmesh,grids.u.inner.ymesh);
	[rhs.inner.u,filtering.u.inner.bcio] = par.bcfunc(grids.u.inner,filtering.u.inner,rhs.inner.u,par,'u');
	filtering.u.inner.bciofull = logical(filtering.u.inner.filterMat'*(1*filtering.u.inner.bcio));
	
	rhs.inner.v = par.rhfunc(grids.v.inner.xmesh,grids.v.inner.ymesh);
	[rhs.inner.v,filtering.v.inner.bcio] = par.bcfunc(grids.v.inner,filtering.v.inner,rhs.inner.v,par,'v');
	filtering.v.inner.bciofull = logical(filtering.v.inner.filterMat'*(1*filtering.v.inner.bcio));
	
	rhs.inner.q = par.rhfunc(grids.q.inner.xmesh,grids.q.inner.ymesh);
	[rhs.inner.q,filtering.q.inner.bcio] = par.bcfunc(grids.q.inner,filtering.q.inner,rhs.inner.q,par,'q');
	filtering.q.inner.bciofull = logical(filtering.q.inner.filterMat'*(1*filtering.q.inner.bcio));
	
	%outer
	rhs.outer.p = par.rhfunc(grids.p.outer.xmesh,grids.p.outer.ymesh);
	[rhs.outer.p,filtering.p.outer.bcio] = par.bcfunc(grids.p.outer,filtering.p.outer,rhs.outer.p,par,'p');
	filtering.p.outer.bciofull = logical(filtering.p.outer.filterMat'*(1*filtering.p.outer.bcio));
	
	rhs.outer.u = par.rhfunc(grids.u.outer.xmesh,grids.u.outer.ymesh);
	[rhs.outer.u,filtering.u.outer.bcio] = par.bcfunc(grids.u.outer,filtering.u.outer,rhs.outer.u,par,'u');
	filtering.u.outer.bciofull = logical(filtering.u.outer.filterMat'*(1*filtering.u.outer.bcio));
	
	rhs.outer.v = par.rhfunc(grids.v.outer.xmesh,grids.v.outer.ymesh);
	[rhs.outer.v,filtering.v.outer.bcio] = par.bcfunc(grids.v.outer,filtering.v.outer,rhs.outer.v,par,'v');
	filtering.v.outer.bciofull = logical(filtering.v.outer.filterMat'*(1*filtering.v.outer.bcio));
	
	rhs.outer.q = par.rhfunc(grids.q.outer.xmesh,grids.q.outer.ymesh);
	[rhs.outer.q,filtering.q.outer.bcio] = par.bcfunc(grids.q.outer,filtering.q.outer,rhs.outer.q,par,'q');
	filtering.q.outer.bciofull = logical(filtering.q.outer.filterMat'*(1*filtering.q.outer.bcio));
	
	if(~exist('figs','var'))
		[grids,filtering,res,par,figs] = par.nssolver(par,grids,filtering,rhs);
	else
		[grids,filtering,res,par,figs] = par.nssolver(par,grids,filtering,rhs,figs);
	end
	
else
	
	[grids,filtering,par] = par.gridmaker(par);
	
	rhs = [];
	
	rhs.inner = par.rhfunc(grids.inner.xmesh,grids.inner.ymesh);
	[rhs.inner,filtering.inner.bcio] = par.bcfunc(grids.inner,filtering.inner,rhs.inner,par);
	filtering.inner.bciofull = logical(filtering.inner.filterMat'*(1*filtering.inner.bcio));
	
	rhs.outer = par.rhfunc(grids.outer.xmesh,grids.outer.ymesh);
	[rhs.outer,filtering.outer.bcio] = par.bcfunc(grids.outer,filtering.outer,rhs.outer,par);
	filtering.outer.bciofull = logical(filtering.outer.filterMat'*(1*filtering.outer.bcio));
	
	[grids,filtering,res,par,figs] = par.nssolver(par,grids,filtering,rhs);
	
end



