function [rhs,filtering] = BCSymChNS(grids,filtering,rhs,par)
	
	[rhs.q.inner ,filtering] = bcq(grids,filtering,rhs.q.inner,par,'inner','inner');
	[rhs.q.outer1,filtering] = bcq(grids,filtering,rhs.q.outer1,par,'outer1','inner');
	[rhs.q.outer2,filtering] = bcq(grids,filtering,rhs.q.outer2,par,'outer2','inner');
	
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
	filtering.q.(side).bciofull = logical(filtering.q.(side).filterMat'*(1*filtering.q.(side).bcio));
	
end

