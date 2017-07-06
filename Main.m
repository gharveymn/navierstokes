if(~exist('par','var'))
	par = Parameters;
end


if(par.usestagger)
	
	[grids,filtering,par] = MakeStaggeredGrids(par);
	
	%implement external force function (on rhs)
	prhs = par.rhfunc(grids.p.outer.xmesh,grids.p.outer.ymesh);
	[prhs,pbcio] = par.bcfunc(grids.p.outer,filtering.p.outer,prhs,par);
	filtering.p.bcio = pbcio;
	rhs.p = prhs;
	
	urhs = par.rhfunc(grids.u.outer.xmesh,grids.u.outer.ymesh);
	[urhs,ubcio] = par.bcfunc(grids.u.outer,filtering.u.outer,urhs,par);
	filtering.u.bcio = ubcio;
	rhs.u = urhs;
	
	vrhs = par.rhfunc(grids.v.outer.xmesh,grids.v.outer.ymesh);
	[vrhs,vbcio] = par.bcfunc(grids.v.outer,filtering.v.outer,vrhs,par);
	filtering.v.bcio = vbcio;
	rhs.v = vrhs;
	
	qrhs = par.rhfunc(grids.q.outer.xmesh,grids.q.outer.ymesh);
	[qrhs,qbcio] = par.bcfunc(grids.q.outer,filtering.q.outer,qrhs,par);
	filtering.q.bcio = qbcio;
	rhs.q = qrhs;
	
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



