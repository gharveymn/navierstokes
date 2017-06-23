function [figs,mat,vec] = InPost(psimesh,bc,grids,filtering,par,figs)
	%INPOST does the post processing of calculation
	
	h = par.h;
	psimeshfull = filtering{1}'*psimesh;
	
	nx = grids{9};
	ny = grids{10};
	
	bcxdfull = bc{1}{2}{1};
	bcydfull = bc{1}{2}{2};
	
	%TODO change so that derivatives are cool at boundaries
	if(par.ghostpoints)

		if(par.filter)
			
			%again we're setting a hard limit of order 3
			switch par.order
				case 1
					%do nothing
				case 2
					[~,~,~,bcxdfull] = closure(grids,filtering,'inner',filtering{5}{1},bcxdfull);
					[~,~,~,bcydfull] = closure(grids,filtering,'inner',filtering{5}{1},bcydfull);
					[~,~,~,psimeshfull] = closure(grids,filtering,'inner',filtering{5}{1},psimeshfull);
					[~,grids,filtering] = closure(grids,filtering,'inner',filtering{5}{1},filtering{5}{1});
				case 3
					[~,~,~,bcxdfull] = closure(grids,filtering,'inner',filtering{5}{2},bcxdfull);
					[~,~,~,bcydfull] = closure(grids,filtering,'inner',filtering{5}{2},bcydfull);
					[~,~,~,psimeshfull] = closure(grids,filtering,'inner',filtering{5}{2},psimeshfull);
					[~,newgrids,newfiltering,gp] = closure(grids,filtering,'inner',filtering{5}{2},filtering{5}{1});

				% 	[~,~,~,umeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
				% 	[~,~,~,vmeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
				% 	[~,~,~,psimeshfull] = closure(newgrids,newfiltering,h,'inner',gp,psimeshfull);
				% 	[~,newgrids,newfiltering,gp] = closure(newgrids,newfiltering,h,'inner',gp,gp);

					[~,~,~,bcxdfull] = closure(newgrids,newfiltering,'inner',gp,bcxdfull);
					[~,~,~,bcydfull] = closure(newgrids,newfiltering,'inner',gp,bcydfull);
					[~,~,~,psimeshfull] = closure(newgrids,newfiltering,'inner',gp,psimeshfull);
					[~,grids,filtering] = closure(newgrids,newfiltering,'inner',gp,gp);
				otherwise
					ME = MException('closure:invalidParameterException','Invalid value for par.order');
					throw(ME)
			end
			
		end
		
% 		%TODO figure out how to get back our psi at the right size
		filterMat = filtering{1};
		on = filtering{3}{1};
		
		bcxd = logical(filterMat*(1*bcxdfull));
		bcyd = logical(filterMat*(1*bcydfull));
		psimesh = filterMat*psimeshfull;
		nx = grids{9};
		ny = grids{10};
		
		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		dx = spdiag(~(bcxd|bcyd))*dx;
		
		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		dy = spdiag(~(bcxd|bcyd))*dy;

		
	else
		filterMat = filtering{1};
		valind = filtering{2};
		on = filtering{3}{1};
		onfull = filtering{3}{2};
		
		%Switch to first order on the boundary
		dbc = boundarysides(grids,filtering);
		bcw = dbc{1};
		bce = dbc{2};
		bcs = dbc{3};
		bcn = dbc{4};
		bcc = dbc{5};

		Dx = sptoeplitz([0 -1],[0 1],nx)./(2*h);
		dx = kron(speye(ny),Dx);
		dx = filterMat*dx*filterMat';
		
		if(par.zeroout)
			dx = spdiag(~(bcw|bce|bcc))*dx;
		else
			dx = spdiag(~bcw)*dx + 1/h*(-spdiag(bcw) + spdiag(bcw(1:end-1),1));
			dx = spdiag(~bce)*dx + 1/h*(-spdiag(bce(2:end),-1) + spdiag(bce));
		end

		Dy = sptoeplitz([0 -1],[0 1],ny)./(2*h);
		dy = kron(Dy,speye(nx));
		dy = filterMat*dy*filterMat';
		
		if(par.zeroout)
			dy = spdiag(~(bcs|bcn|bcc))*dy;
		else
			dy = spdiag(~bcs)*dy + 1/h*(-spdiag(bcs) + spdiag(bcs(1:end-nx),nx));
			dy = spdiag(~bcn)*dy + 1/h*(-spdiag(bcn(nx+1:end),-nx) + spdiag(bcn));
		end
	end
	
	umesh = dy*psimesh;
	vmesh = -dx*psimesh;
	
	umeshfull = filterMat'*umesh;
	Umesh = reshape(umeshfull,[nx,ny])';
	
	vmeshfull = filterMat'*vmesh;
	Vmesh = reshape(vmeshfull,[nx,ny])';
	
	psimeshfull = filterMat'*psimesh;
	Psimesh = reshape(psimeshfull,[nx,ny])';
	
	% use matrices rather than cell arrays so they throw a dimension error if we have a bug mismatched
	mat = cat(3,grids{5},grids{6},Umesh,Vmesh,Psimesh);
	vec = cat(2,grids{3},grids{4},umesh,vmesh,psimesh);
	
	if(par.plot)
		if(exist('figs','var'))
			figs = Plot(mat,vec,par,figs);
		else
			figs = Plot(mat,vec,par);
		end
	end
	
end

