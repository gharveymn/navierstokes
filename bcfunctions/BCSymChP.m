function [rhs,bc] = BCSymChP(grids,filtering,rhs,par)
	%BCSYMCHP BCSymCh for the primitive formulation
	
	xmesh = grids{3};
	ymesh = grids{4};
	nx = grids{9};
	valindinner = filtering{2}{1};
	valindouter = filtering{2}{2};
	
	on = filtering{3}{1};
	bcfull = filtering{4}{2};
	gpca = filtering{5};
	
	del = par.h;
	
	%for use with symch map
	bc = {on,on,on};
	
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	ymax = max(ymesh(on));
	ymin = min(ymesh(on));
	
	
	% inflow
	inflowmax = max(ymesh(xmesh==xmin & on));
	inflowmin = min(ymesh(xmesh==xmin & on));
	
	h = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + h;
	
	a = 1;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(h^2 - (ymesh-centerin).^2);
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs{1} = rhs{1} + in;
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax & on));
	outflowmin = min(ymesh(xmesh==xmax & on));
	
	H = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + H;
	
	e = (4*a*h^3/3)*3/(4*H^3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(H^2 - (ymesh-centerout).^2);
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs{1} = rhs{1} + out;
	
	% set top
	rhs{1}(ymesh > centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymax);
	
	% set bottom
	rhs{1}(ymesh < centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymin);
	
		if(par.order > 1)
		rhs = extendgp(rhs,dbcfull,valindouter,gpca,nx);
	end
	
	for i=1:par.order-1
		bc{1}{1} = bc{1}{1}|gpca{i}(valindouter);
		bc{1}{2} = bc{1}{1};
		bc{2}{1} = logical(filtering{1}'*(1*bc{1}{1}));
		bc{2}{2} = logical(filtering{1}'*(1*bc{1}{2}));
	end
	
	if(par.ghostpoints)
		rhs{1} = extendgp(rhs{1},bcfull,valindouter,gpca,nx);
		bc{1} = bc{1}|gpca{1}(valindouter)|gpca{2}(valindouter);	
		bc{3} = bc{1};
	end
	bc{1} = bc{1}&((xmesh<=inflowx & ymesh <= inflowmax & ymesh >= inflowmin)...
						|(xmesh>=outflowx & ymesh <= outflowmax & ymesh >= outflowmin));				
	
	if(par.ghostpoints)
		bc{2} = (~bc{1}&(gpca{1}(valindouter)|gpca{2}(valindouter)|on))...
								|((xmesh==outflowx|xmesh==outflowx+del|xmesh==outflowx+2*del)&ymesh==outflowmin)...
								|((xmesh==outflowx|xmesh==outflowx+del|xmesh==outflowx+2*del)&ymesh==outflowmax)...
								|((xmesh==inflowx|xmesh==inflowx-del|xmesh==inflowx-2*del)&ymesh==inflowmin)...
								|((xmesh==inflowx|xmesh==inflowx-del|xmesh==inflowx-2*del)&ymesh==inflowmax);
	else
		bc{2} = (~bc{1}&on)|(xmesh==outflowx&ymesh==outflowmin)...
								|(xmesh==outflowx&ymesh==outflowmax)...
								|(xmesh==inflowx&ymesh==inflowmin)...
								|(xmesh==inflowx&ymesh==inflowmax);
	end
	
end

