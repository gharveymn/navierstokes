if(~exist('par','var'))
	par = Parameters;
end

if(~exist('figs','var'))
	[figs,mat,vec] = run(par);
else
	[figs,mat,vec] = run(par,figs);
end

X = mat(:,:,1);
Y = mat(:,:,2);
U = mat(:,:,3);
V = mat(:,:,4);
Psi = mat(:,:,5);

x = vec(:,1);
y = vec(:,2);
u = vec(:,3);
v = vec(:,4);
psi = vec(:,5);

clear mat vec

if(~exist('testrun','var') || ~testrun)
	clear par
end

function [figs,mat,vec] = run(par,figs)
	
	rhfunc = par.rhfunc;
	bcfunc = par.bcfunc;
	solver = par.solver;
	ddsolver = par.ddsolver;
	
	[grids,filtering,par] = ParseValidIndices(par);
	
	nx = numel(grids{1});
	ny = numel(grids{2});
	
	%implement external force function (on rhs)
	rhs = rhfunc(grids{3},grids{4});
	
	%make right hand side for Dirichlet BCs and get indices for those points
	[rhs,bc] = bcfunc(grids,filtering,rhs,par);
	
	%  	filterMat = filtering{1};
	%   	rmeshfull = filterMat'*rhs;
	%   	Rmesh = reshape(rmeshfull,[nx,ny])';
	%   	surf(grids{5},grids{6},Rmesh,'edgecolor','none','facecolor','interp');
	% 	clr = abs(rmeshfull)./norm(rmeshfull(isfinite(rmeshfull)),inf);
	% 	clrs = [clr zeros(numel(clr),1) 1-clr];
	%   	scatter3(grids{7},grids{8},rmeshfull,[],clrs,'.');
	
	filterMat = filtering{1};
	
	if(par.streamfunction)
		if(par.ddrun)
			if(exist('figs','var'))
				psimesh = ddsolver(grids,filtering,rhs,bc,par,solver,figs);
			else
				psimesh = ddsolver(grids,filtering,rhs,bc,par,solver);
			end
		else
			psimesh = solver(grids,filtering,rhs,bc);
		end
		
		if(exist('figs','var'))
			[figs,mat,vec] = InPost(psimesh,bc,grids,filtering,par,figs);
		else
			[figs,mat,vec] = InPost(psimesh,bc,grids,filtering,par);
		end
	else
		
		[umesh,vmesh,pmesh] = solver(nx,ny,bc,rhs,filterMat,h);
		
		umeshfull = filterMat'*umesh;
		Umesh = reshape(umeshfull,[nx,ny])';
		
		vmeshfull = filterMat'*vmesh;
		Vmesh = reshape(vmeshfull,[nx,ny])';
		
		psimeshfull = filterMat'*pmesh;
		Pmesh = reshape(psimeshfull,[nx,ny])';
		
		if(par.filter)
			on = filtering{3}{1};
			grids{3} = grids{3}(~on);
			grids{4} = grids{4}(~on);
			umesh = umesh(~on);
			vmesh = vmesh(~on);
			pmesh = pmesh(~on);
		end
		
		mat = cat(3,grids{5},grids{6},Umesh,Vmesh,Pmesh);
		vec = cat(2,grids{3},grids{4},umesh,vmesh,pmesh);
		
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

