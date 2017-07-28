function [rhs,bcio] = BCSymChNS(grids,filtering,rhs,par,grid)
	switch grid
		case 'u'
			[rhs,bcio] = bcu(grids,filtering,rhs,par);
		case 'v'
			[rhs,bcio] = bcv(grids,filtering,rhs,par);
		case 'p'
			[rhs,bcio] = bcp(grids,filtering,rhs,par);
		case 'q'
			[rhs,bcio] = bcq(grids,filtering,rhs,par);
		otherwise
			ME = MException('bcsymchns:invalidParameterException','Invalid value for grid');
			throw(ME)
	end
end

function [rhs,bcio] = bcu(grids,filtering,rhs,par)
	
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	
	on = filtering.on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	% inflow
	inflowmax = max(ymesh(xmesh==xmin & on));
	inflowmin = min(ymesh(xmesh==xmin & on));
	
	d = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + d;
	
	a = par.inflowAmp;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(d^2 - (ymesh-centerin).^2);
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax & on));
	outflowmin = min(ymesh(xmesh==xmax & on));
	
	D = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + D;
	
	e = (4*a*d^3/3)*3/(4*D^3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(D^2 - (ymesh-centerout).^2);
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(ymesh > centerout & ~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = 0;
	
	% set bottom
	rhs(ymesh < centerout & ~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = 0;
	
	bcio = (xmesh==inflowx|xmesh==outflowx)&on;
end

function [rhs,bcio] = bcv(grids,filtering,rhs,par)
	
	xmesh = grids.xmesh;	
	on = filtering.on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	
end

function [rhs,bcio] = bcp(grids,filtering,rhs,par)
	
	xmesh = grids.xmesh;	
	on = filtering.on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	
end

function [rhs,bcio] = bcq(grids,filtering,rhs,par)
	
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	
	on = filtering.on;	
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	% inflow
	inflowmax = max(ymesh(xmesh==xmin & on));
	inflowmin = min(ymesh(xmesh==xmin & on));
	
	d = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + d;
	
	a = par.inflowAmp;
	c = a*(d^2.*(inflowmax-centerin) - (inflowmax-centerin).^3./3);
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(d^2.*(ymesh-centerin) - (ymesh-centerin).^3./3) + c;
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax & on));
	outflowmin = min(ymesh(xmesh==xmax & on));
	
	D = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + D;
	
	e = (4*a*d^3/3)*3/(4*D^3);
	f = e*(D^2.*(outflowmax-centerout) - (outflowmax-centerout).^3./3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(D^2.*(ymesh-centerout) - (ymesh-centerout).^3./3) + f;
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(ymesh > centerout & ~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = max(rhs(xmesh==inflowx));
	
	% set bottom
	rhs(ymesh < centerout & ~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = min(rhs(xmesh==inflowx));
	
	bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	
end

