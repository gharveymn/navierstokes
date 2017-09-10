function [rhs,filtering] = BCConIO(grids,filtering,rhs,par)
	
	[rhs.q.inner ,filtering] = bcq(grids,filtering,rhs.q.inner ,par,'inner','inner');
	[rhs.q.outer1,filtering] = bcq(grids,filtering,rhs.q.outer1,par,'outer1','inner');
	[rhs.q.outer2,filtering] = bcq(grids,filtering,rhs.q.outer2,par,'outer2','inner');
	
end

function [rhs,filtering] = bcq(grids,filtering,rhs,par,side)
	
	xmesh = grids.q.(side).xmesh;
	
	on = filtering.q.(side).on;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	
	inflowx = xmin*ones(numel(xmesh),1);
	outflowx = xmax*ones(numel(xmesh),1);
	
	filtering.q.(side).bcio = ((xmesh==inflowx) & on) | ((xmesh==outflowx) & on);
	filtering.q.(side).bciofull = logical(filtering.q.(side).filterMat'*(1*filtering.q.(side).bcio));
	
end

