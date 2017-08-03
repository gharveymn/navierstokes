function rhs = RHZero(grids)

	rhs = [];

	rhs.u.inner = zeros(numel(grids.u.inner.xmesh),1);
	rhs.v.inner = zeros(numel(grids.v.inner.xmesh),1);
	rhs.p.inner = zeros(numel(grids.p.inner.xmesh),1);
	rhs.q.inner = zeros(numel(grids.q.inner.xmesh),1);
	
	rhs.u.outer = zeros(numel(grids.u.outer.xmesh),1);
	rhs.v.outer = zeros(numel(grids.v.outer.xmesh),1);
	rhs.p.outer = zeros(numel(grids.p.outer.xmesh),1);
	rhs.q.outer = zeros(numel(grids.q.outer.xmesh),1);
end

