function rhs = RHOne(grids,par)
	
	for i=1:numel(par.varnames)
		
		rhs.(par.varnames{i}).inner = ones(numel(grids.(par.varnames{i}).inner.xmesh),1);
		rhs.(par.varnames{i}).outer = ones(numel(grids.(par.varnames{i}).outer.xmesh),1);
	
	end
	
end

