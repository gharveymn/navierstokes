if(~exist('testrun','var') || ~testrun)
	par = Parameters;
end


if(par.usestagger)
	
	[grids,filtering,par] = par.gridmaker(par);
	
	rhs = par.rhfunc(grids);
	[rhs,filtering] = par.bcfunc(grids,filtering,rhs,par);
	
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



