function [rhs,filtering] = BCSymChNS(grids,filtering,rhs,par,spec)
	if(~exist('spec','var'))
		for i=1:numel(par.varnames)
			
			bcf = str2func(['bc' par.varnames{i}]);
			
			[rhs.(par.varnames{i}).inner,filtering] = bcf(grids,filtering,rhs.(par.varnames{i}).inner,par,'inner');
			[rhs.(par.varnames{i}).outer,filtering] = bcf(grids,filtering,rhs.(par.varnames{i}).outer,par,'outer');
			
		end
	else
		bcf = str2func(['bc' spec.varname]);
		[rhs.(spec.varname).(spec.side),filtering] = bcf(grids,filtering,rhs.(spec.varname).(spec.side),par,spec.side);
	end
	
end



function [rhs,filtering] = bcu(grids,filtering,rhs,par,side)
	
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
	
	filtering.u.(side).bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	filtering = makebciofull(filtering,'u',side);
end

function [rhs,filtering] = bcv(grids,filtering,rhs,par,side)
	
	xmesh = grids.v.(side).xmesh;
	on = filtering.v.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	filtering.v.(side).bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	filtering = makebciofull(filtering,'v',side);
	
end

function [rhs,filtering] = bcp(grids,filtering,rhs,par,side)
	
	xmesh = grids.p.(side).xmesh;
	on = filtering.p.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	filtering.p.(side).bcio = (xmesh==inflowx|xmesh==outflowx)&on;
	filtering = makebciofull(filtering,'p',side);
	
end

function [rhs,filtering] = bcq(grids,filtering,rhs,par,side,sideinner)
	
	xmesh = grids.q.(side).xmesh;
	ymesh = grids.q.(side).ymesh;
	
	on = filtering.q.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	% inflow
	if(exist('sideinner','var'))
		ixmax = max(grids.q.(sideinner).xmesh(filtering.q.(sideinner).on));
		ixmin = min(grids.q.(sideinner).xmesh(filtering.q.(sideinner).on));
		
		inflowmax = max(grids.q.(sideinner).ymesh(grids.q.(sideinner).xmesh==ixmin & filtering.q.(sideinner).on));
		inflowmin = min(grids.q.(sideinner).ymesh(grids.q.(sideinner).xmesh==ixmin & filtering.q.(sideinner).on));
	else
		inflowmax = max(ymesh(xmesh==xmin & on));
		inflowmin = min(ymesh(xmesh==xmin & on));
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
	if(exist('sideinner','var'))
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
	rhsmax = max(rhs(xmesh==inflowx));
	rhs(ymesh > centerout & (~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) | ((ymesh>inflowmax)&~(xmesh==outflowx)) | ((ymesh>outflowmax)&~(xmesh==inflowx))) & on) = rhsmax;
	
	% set bottom
	rhsmin = min(rhs(xmesh==inflowx));
	rhs(ymesh < centerout & (~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx)) | ((ymesh<inflowmin)&~(xmesh==outflowx)) | ((ymesh<outflowmin)&~(xmesh==inflowx))) & on) = rhsmin;
	
	filtering.q.(side).bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	filtering = makebciofull(filtering,'q',side);
	
end

function [filtering] = makebciofull(filtering,varname,side)
	filtering.(varname).(side).bciofull = logical(filtering.(varname).(side).filterMat'*(1*filtering.(varname).(side).bcio));
end

