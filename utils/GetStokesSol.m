function [figs,mat,vec,grids,filtering,bc,rhs] = GetStokesSol(par,figs)
	
	rhfunc = par.rhfunc;
	bcfunc = par.bcfunc;
	solver = par.solver;
	ddsolver = par.ddsolver;
	
	[grids,filtering,par] = MakeGrids(par);
	
	nxp1 = numel(grids.xinit);
	nyp1 = numel(grids.yinit);
	
	%implement external force function (on rhs)
	rhs = rhfunc(grids.xmesh,grids.ymesh);
	
	%make right hand side for Dirichlet BCs and get indices for those points
	[rhs,bc] = bcfunc(grids,filtering,rhs,par);
	
	%  	filterMat = filtering.filterMat;
	%   	rmeshfull = filterMat'*rhs;
	%   	Rmesh = reshape(rmeshfull,[nx,ny])';
	%   	surf(grids.Xmesh,grids.Ymesh,Rmesh,'edgecolor','none','facecolor','interp');
	% 	clr = abs(rmeshfull)./norm(rmeshfull(isfinite(rmeshfull)),inf);
	% 	clrs = [clr zeros(numel(clr),1) 1-clr];
	%   	scatter3(grids.xmeshfull,grids.ymeshfull,rmeshfull,[],clrs,'.');
	
	filterMat = filtering.filterMat;
	
	if(par.streamfunction)
		if(par.ddrun)
			if(exist('figs','var'))
				qmesh = ddsolver(grids,filtering,rhs,bc,par,solver,figs);
			else
				qmesh = ddsolver(grids,filtering,rhs,bc,par,solver);
			end
		else
			qmesh = solver(grids,filtering,rhs,bc);
		end
		
		if(exist('figs','var'))
			[figs,mat,vec] = InPost(qmesh,bc,grids,filtering,par,figs);
		else
			[figs,mat,vec] = InPost(qmesh,bc,grids,filtering,par);
		end
	else
		
		[umesh,vmesh,pmesh] = solver(nxp1,nyp1,bc,rhs,filterMat,h);
		
		umeshfull = filterMat'*umesh;
		Umesh = reshape(umeshfull,[nxp1,nyp1])';
		
		vmeshfull = filterMat'*vmesh;
		Vmesh = reshape(vmeshfull,[nxp1,nyp1])';
		
		qmeshfull = filterMat'*pmesh;
		Pmesh = reshape(qmeshfull,[nxp1,nyp1])';
		
		if(par.filter)
			on = filtering.on;
			grids.xmesh = grids.xmesh(~on);
			grids.ymesh = grids.ymesh(~on);
			umesh = umesh(~on);
			vmesh = vmesh(~on);
			pmesh = pmesh(~on);
		end
		
		mat = cat(3,grids.Xmesh,grids.Ymesh,Umesh,Vmesh,Pmesh);
		vec = cat(2,grids.xmesh,grids.ymesh,umesh,vmesh,pmesh);
		
		if(par.plot)
			if(exist('figs','var'))
				figs = Plot(mat,vec,par,figs);
			else
				figs = Plot(mat,vec,par);
			end
		else
			figs = [];
		end
		
	end
	
	
end