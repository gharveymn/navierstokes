function [rhs,bcinds] = BCSymCh(grids,filtering,rhs,par)
	
	xmeshfull = grids.inner.xmeshfull;
	ymeshfull = grids.inner.ymeshfull;
	xmesh = grids.inner.xmesh;
	ymesh = grids.inner.ymesh;
	nx = grids.nxp1;
	valindinner = filtering.valindinner;
	valindouter = filtering.valindouter;
	
	on = filtering.on;
	bcfull = filtering.dbcfull;
	gpca = filtering.gp;
	
	del = par.h;
	
	%for use with symch map
	bcinds = 0*xmesh;
	
	% add all the the indices which are on the boundary
	bcinds = bcinds | on;
	
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
	c = 1/12;
	
	inflowx = xmin*ones(numel(xmesh),1);
	in = a*(h^2.*(ymesh-centerin) - (ymesh-centerin).^3./3) + c;
	in(~(xmesh==inflowx) | ~on) = 0;
	
	rhs = rhs + in;
	
	%outflow
	outflowmax = max(ymesh(xmesh==xmax & on));
	outflowmin = min(ymesh(xmesh==xmax & on));
	
	H = (outflowmax-outflowmin)/2;
	
	centerout = outflowmin + H;
	
	m = (max(in)-min(in))/(ymax-ymin);
	
	outflowx = xmax*ones(numel(xmesh),1);
	out = m*(ymesh-ymin) + min(in);
	out(~(xmesh==outflowx) | ~on) = 0;
	
	rhs = rhs + out;
	
	% set top
	rhs(ymesh > centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymax);
	
	% set bottom
	rhs(ymesh < centerout & ~(xmesh==inflowx | xmesh==outflowx) & on) = out(xmesh==xmax&ymesh==ymin);
	
	rhs = extendgp(rhs,bcfull,valindouter,gpca,nx);
	bcinds = bcinds|gpca{1}(valindouter)|gpca{2}(valindouter);
	
end

