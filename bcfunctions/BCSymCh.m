function [rhs,bcio] = BCSymCh(grids,filtering,rhs,par)
	
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	nxp1 = par.nxp1;
	valind = filtering.valind;
	
	on = filtering.on;
	onfull = filtering.onfull;
	dbcfull = filtering.dbcfull;
	gpca = filtering.gp;
	
	h = par.h;
	
	% add all the the indices which are on the boundary
	% bc{1} - dirichlet
	% bc{2} - neumann
	bcio = {{{on,on},{onfull,onfull}},{{[],[]},{[],[]}}};
	
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
	rhs(ymesh > centerout & ~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = 0;
	
	% set bottom
	rhs(ymesh < centerout & ~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = 0;
	
	if(par.ghostpoints)
		if(par.order > 1)
			rhs = extendgp(rhs,dbcfull,valind,gpca,nxp1);
		end

		for i=1:par.order-1
			bcio{1}{1}{1} = bcio{1}{1}{1}|gpca{i}(valind);
			bcio{1}{1}{2} = bcio{1}{1}{1};
			bcio{1}{2}{1} = logical(filtering.filterMat'*(1*bcio{1}{1}{1}));
			bcio{1}{2}{2} = logical(filtering.filterMat'*(1*bcio{1}{1}{2}));
		end
	else
		bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	end
	
	filtering.bcinout = bcio;
	
end

