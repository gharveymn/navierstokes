function rhs = RHZero(grids)

	rhs.q.inner = zeros(numel(grids.q.inner.xmesh),1);
	rhs.q.outer1 = zeros(numel(grids.q.outer1.xmesh),1);
	rhs.q.outer2 = zeros(numel(grids.q.outer2.xmesh),1);
	
end

