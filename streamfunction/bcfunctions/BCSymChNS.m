function [rhs,filtering] = BCSymChNS(grids,filtering,rhs,par)
	
	if(par.stagger)
		[rhs,filtering] = getBC(grids,filtering,rhs,par);
	else
		[rhs.q.inner,filtering.q.inner.bcio] = bcq(grids,filtering,rhs.q.inner,par,'inner','inner');
		filtering.q.inner.bciofull = logical(filtering.q.inner.filterMat'*(1*filtering.q.inner.bcio));
		
		[rhs.q.outer1,filtering.q.outer1.bcio] = bcq(grids,filtering,rhs.q.outer1,par,'outer1','inner');
		filtering.q.outer1.bciofull = logical(filtering.q.outer1.filterMat'*(1*filtering.q.outer1.bcio));
		
		[rhs.q.outer2,filtering.q.outer2.bcio] = bcq(grids,filtering,rhs.q.outer2,par,'outer2','inner');
		filtering.q.outer2.bciofull = logical(filtering.q.outer2.filterMat'*(1*filtering.q.outer2.bcio));
	end
	
end

function [rhs,filtering] = getBC(grids,filtering,rhs,par,spec)
	if(~exist('spec','var'))
		for i=1:numel(par.varnames)
			
			bcf = str2func(['bc' par.varnames{i}]);
			
			[rhs.(par.varnames{i}).inner,filtering.(par.varnames{i}).inner.bcio] = bcf(grids,filtering,rhs.(par.varnames{i}).inner,par,'inner');
			filtering.(par.varnames{i}).inner.bciofull = logical(filtering.(par.varnames{i}).inner.filterMat'*(1*filtering.(par.varnames{i}).inner.bcio));
			
			[rhs.(par.varnames{i}).outer,filtering.(par.varnames{i}).outer.bcio] = bcf(grids,filtering,rhs.(par.varnames{i}).outer,par,'outer');
			filtering.(par.varnames{i}).outer.bciofull = logical(filtering.(par.varnames{i}).outer.filterMat'*(1*filtering.(par.varnames{i}).outer.bcio));
			
		end
	else
		bcf = str2func(['bc' spec.varname]);
		
		[rhs.(spec.varname).(spec.side),filtering.(spec.varname).(spec.side).bcio] = bcf(grids,filtering,rhs.(spec.varname).(spec.side),par,spec.side);
		filtering.(spec.varname).(spec.side).bciofull = logical(filtering.(spec.varname).(spec.side).filterMat'*(1*filtering.(spec.varname).(spec.side).bcio));
	end
end



function [rhs,bcio] = bcu(grids,filtering,rhs,par,side)
	
	xmesh = grids.u.(side).xmesh;
	ymesh = grids.u.(side).ymesh;
	on = filtering.u.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	% inflow
	inflowmax = max(grids.u.inner.ymesh(grids.u.inner.xmesh==min(grids.u.inner.xmesh(filtering.u.inner.on)) & filtering.u.inner.on));
	inflowmin = min(grids.u.inner.ymesh(grids.u.inner.xmesh==min(grids.u.inner.xmesh(filtering.u.inner.on)) & filtering.u.inner.on));
	
	d = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + d;
	
	a = par.inflowAmp;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(d^2 - (ymesh-centerin).^2);
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(grids.u.inner.ymesh(grids.u.inner.xmesh==max(grids.u.inner.xmesh(filtering.u.inner.on)) & filtering.u.inner.on));
	outflowmin = min(grids.u.inner.ymesh(grids.u.inner.xmesh==max(grids.u.inner.xmesh(filtering.u.inner.on)) & filtering.u.inner.on));
	
	D = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + D;
	
	e = (4*a*d^3/3)*3/(4*D^3);
	f = e*(D^2.*(outflowmax-centerout) - (outflowmax-centerout).^3./3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = (e*(D^2.*(outflowmax-centerout) - (outflowmax-centerout).^3./3) + f)/(2*D)*ones(numel(xmesh),1);
	out(~(xmesh==outflowx) | ~on) = 0;
	rhs = rhs + out;
	
	% set rest
	if(strcmp(side,'outer'))
		rhs((~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx))|filtering.u.(side).dbc.c) & on) = 0;
	else
		rhs(~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) & on) = 0;
	end
	
	bcio = (xmesh==inflowx|xmesh==outflowx)&on;
end

function [rhs,bcio] = bcv(grids,filtering,rhs,par,side)
	
	xmesh = grids.v.(side).xmesh;
	on = filtering.v.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	
end

function [rhs,bcio] = bcp(grids,filtering,rhs,par,side)
	
	xmesh = grids.p.(side).xmesh;
	on = filtering.p.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	
end

function [rhs,bcio] = bcq(grids,filtering,rhs,par,side,sideinner)
	
	xmesh = grids.q.(side).xmesh;
	ymesh = grids.q.(side).ymesh;
	
	on = filtering.q.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	% inflow
	if(exist('sideouter','var'))
		ixmax = max(grids.q.(sideinner).xmesh(filtering.q.(sideinner).on));
		ixmin = min(grids.q.(sideinner).xmesh(filtering.q.(sideinner).on));
		
		inflowmax = max(grids.q.(sideinner).ymesh(grids.q.(sideinner).xmesh==ixmin & filtering.q.(sideinner).on));
		inflowmin = min(grids.q.(sideinner).ymesh(grids.q.(sideinner).xmesh==ixmin & filtering.q.(sideinner).on));
	else
		inflowmax = max(grids.q.(sideinner).ymesh(xmesh==xmin & on));
		inflowmin = min(grids.q.(sideinner).ymesh(xmesh==xmin & on));
	end
	
	d = (inflowmax-inflowmin)/2;
	
	centerin = inflowmin + d;
	
	a = par.inflowAmp;
	c = a*(d^2.*(inflowmax-centerin) - (inflowmax-centerin).^3./3);
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(d^2.*(ymesh-centerin) - (ymesh-centerin).^3./3) + c;
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	if(exist('sideouter','var'))
		outflowmax = max(grids.q.(sideinner).ymesh(grids.q.(sideinner).xmesh==ixmax & filtering.q.(sideinner).on));
		outflowmin = min(grids.q.(sideinner).ymesh(grids.q.(sideinner).xmesh==ixmax & filtering.q.(sideinner).on));
	else
		outflowmax = max(ymesh(xmesh==xmax & on));
		outflowmin = min(ymesh(xmesh==xmax & on));
	end
	
	D = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + D;
	
	e = (4*a*d^3/3)*3/(4*D^3);
	f = e*(D^2.*(outflowmax-centerout) - (outflowmax-centerout).^3./3);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*(D^2.*(ymesh-centerout) - (ymesh-centerout).^3./3) + f;
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(ymesh > centerout & (~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) | ymesh > inflowmax) & on) = max(rhs(xmesh==inflowx));
	
	% set bottom
	rhs(ymesh < centerout & (~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) | ymesh < inflowmin) & on) = min(rhs(xmesh==inflowx));
	
	bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	
end

