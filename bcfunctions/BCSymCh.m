function [rhs,bc] = BCSymCh(grids,filtering,rhs,par)
	
	xmeshfull = grids.xmeshfull;
	ymeshfull = grids.ymeshfull;
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	nx = grids.nx;
	valindinner = filtering.valindinner;
	valindouter = filtering.valindouter;
	
	on = filtering.on;
	onfull = filtering.onfull;
	dbcfull = filtering.dbcfull;
	gpca = filtering.gp;
	
	h = par.h;
	
	% add all the the indices which are on the boundary
	% bc{1} - dirichlet
	% bc{2} - neumann
	bc = {{{on,on},{onfull,onfull}},{{[],[]},{[],[]}}};
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	ymax = max(ymesh(on));
	ymin = min(ymesh(on));
	
	% inflow
	inflowmax = max(ymesh(xmesh==xmin & on));
	inflowmin = min(ymesh(xmesh==xmin & on));
	
	d = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + d;
	
	a = par.inflowAmp;
	c = 1/12;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(d^2.*(ymesh-centerin) - (ymesh-centerin).^3./3) + a*c;
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax & on));
	outflowmin = min(ymesh(xmesh==xmax & on));
	
	D = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + D;
	
	f = 1/12;
	e = (4*a*d^3/3 + c - f)*3/(4*D^3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(D^2.*(ymesh-centerout) - (ymesh-centerout).^3./3) + a*f;
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(ymesh > centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymax);
	
	% set bottom
	rhs(ymesh < centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymin);
	
	if(par.order > 1)
		rhs = extendgp(rhs,dbcfull,valindouter,gpca,nx);
	end
	
	for i=1:par.order-1
		bc{1}{1}{1} = bc{1}{1}{1}|gpca{i}(valindouter);
		bc{1}{1}{2} = bc{1}{1}{1};
		bc{1}{2}{1} = logical(filtering.filterMat'*(1*bc{1}{1}{1}));
		bc{1}{2}{2} = logical(filtering.filterMat'*(1*bc{1}{1}{2}));
	end
	
end

