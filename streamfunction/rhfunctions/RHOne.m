function rhs = RHOne(grids)
	rhs.q.inner = ones(numel(grids.q.inner.xmesh),1);
	rhs.q.outer1 = ones(numel(grids.q.outer1.xmesh),1);
	rhs.q.outer2 = ones(numel(grids.q.outer2.xmesh),1);
end

