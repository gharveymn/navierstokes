function [rhs,filtering] = BCSymChNS(grids,filtering,rhs,par)
	
	for i=1:numel(par.varnames)
		
		bcf = str2func(['bc' par.varnames{i}]);
		
		[rhs.(par.varnames{i}).inner,filtering.(par.varnames{i}).inner.bcio] = bcf(grids,filtering,rhs.(par.varnames{i}).inner,par,'inner');
		filtering.(par.varnames{i}).inner.bciofull = logical(filtering.(par.varnames{i}).inner.filterMat'*(1*filtering.(par.varnames{i}).inner.bcio));
		
		[rhs.(par.varnames{i}).outer,filtering.(par.varnames{i}).outer.bcio] = bcf(grids,filtering,rhs.(par.varnames{i}).outer,par,'outer');
		filtering.(par.varnames{i}).outer.bciofull = logical(filtering.(par.varnames{i}).outer.filterMat'*(1*filtering.(par.varnames{i}).outer.bcio));
	
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
	rhs((~(effeq(xmesh,inflowx) | effeq(xmesh,outflowx))|filtering.u.(side).dbc.c) & on) = 0;
	
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

function [rhs,bcio] = bcq(grids,filtering,rhs,par,side)
	
	xmesh = grids.q.(side).xmesh;
	
	on = filtering.q.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	
end

