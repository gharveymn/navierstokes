function [rhs,filtering] = BCConIO(grids,filtering,rhs,par,spec)
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
	
	a = par.inflowAmp;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*ones(size(xmesh));
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(grids.u.inner.ymesh(grids.u.inner.xmesh==max(grids.u.inner.xmesh(filtering.u.inner.on)) & filtering.u.inner.on));
	outflowmin = min(grids.u.inner.ymesh(grids.u.inner.xmesh==max(grids.u.inner.xmesh(filtering.u.inner.on)) & filtering.u.inner.on));
	
	D = (outflowmax-outflowmin)/2;
	
	e = a*(d/D);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = e*ones(size(xmesh));
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

function [rhs,filtering] = bcq(grids,filtering,rhs,par,side)
	
	xmesh = grids.q.(side).xmesh;
	
	on = filtering.q.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	filtering.q.(side).bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	filtering = makebciofull(filtering,'q',side);
	
end

function [filtering] = makebciofull(filtering,varname,side)
	filtering.(varname).(side).dbc.io = filtering.(varname).(side).bcio;
	filtering.(varname).(side).bciofull = logical(filtering.(varname).(side).filterMat'*(1*filtering.(varname).(side).bcio));
	filtering.(varname).(side).dbcfull.io = filtering.(varname).(side).bciofull;
end

