if(~exist('par','var'))
	par = Parameters;
end


if(par.usestagger)
	
	[grids,filtering,par] = MakeStaggeredGrids(par);
	
	%inner
	rhs.inner.p = par.rhfunc(grids.p.inner.xmesh,grids.p.inner.ymesh);
	[rhs.inner.p,filtering.p.inner.bcio] = par.bcfunc(grids.p.inner,filtering.p.inner,rhs.inner.p,par);
	
	rhs.inner.u = par.rhfunc(grids.u.inner.xmesh,grids.u.inner.ymesh);
	[rhs.inner.u,filtering.u.inner.bcio] = par.bcfunc(grids.u.inner,filtering.u.inner,rhs.inner.u,par);
	
	rhs.inner.v = par.rhfunc(grids.v.inner.xmesh,grids.v.inner.ymesh);
	[rhs.inner.v,filtering.v.inner.bcio] = par.bcfunc(grids.v.inner,filtering.v.inner,rhs.inner.v,par);
	
	rhs.inner.q = par.rhfunc(grids.q.inner.xmesh,grids.q.inner.ymesh);
	[rhs.inner.q,filtering.q.inner.bcio] = par.bcfunc(grids.q.inner,filtering.q.inner,rhs.inner.q,par);
	
	%outer
	rhs.outer.p = par.rhfunc(grids.p.outer.xmesh,grids.p.outer.ymesh);
	[rhs.outer.p,filtering.p.outer.bcio] = par.bcfunc(grids.p.outer,filtering.p.outer,rhs.outer.p,par);
	
	rhs.outer.u = par.rhfunc(grids.u.outer.xmesh,grids.u.outer.ymesh);
	[rhs.outer.u,filtering.u.outer.bcio] = par.bcfunc(grids.u.outer,filtering.u.outer,rhs.outer.u,par);
	
	rhs.outer.v = par.rhfunc(grids.v.outer.xmesh,grids.v.outer.ymesh);
	[rhs.outer.v,filtering.v.outer.bcio] = par.bcfunc(grids.v.outer,filtering.v.outer,rhs.outer.v,par);
	
	rhs.outer.q = par.rhfunc(grids.q.outer.xmesh,grids.q.outer.ymesh);
	[rhs.outer.q,filtering.q.outer.bcio] = par.bcfunc(grids.q.outer,filtering.q.outer,rhs.outer.q,par);
	
	if(~exist('figs','var'))
		[par,mat,vec] = par.nssolver(par,grids,filtering,rhs);
	else
		[par,mat,vec] = par.nssolver(par,grids,filtering,rhs,figs);
	end
	
else
	if(~exist('figs','var'))
		[figs,mat,vec,grids,filtering,bc,rhs,par] = GetStokesSol(par);
	else
		[figs,mat,vec,grids,filtering,bc,rhs,par] = GetStokesSol(par,figs);
	end
	
	[par,mat,vec] = par.nssolver(par,figs,mat,vec,grids,filtering,bc,rhs);
	
end

X = mat(:,:,1);
Y = mat(:,:,2);
U = mat(:,:,3);
V = mat(:,:,4);
Q = mat(:,:,5);

x = vec(:,1);
y = vec(:,2);
u = vec(:,3);
v = vec(:,4);
q = vec(:,5);

clear mat vec

if(~exist('testrun','var') || ~testrun)
	clear par
end



