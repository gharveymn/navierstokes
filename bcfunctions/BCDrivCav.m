function [rhs,filtering] = BCDrivCav(grids,filtering,rhs,par)
	
	for i=1:numel(par.varnames)
		
		bcf = str2func(['bc' par.varnames{i}]);
		
		[rhs.(par.varnames{i}).inner,filtering.(par.varnames{i}).inner.bcio] = bcf(grids,filtering,rhs.(par.varnames{i}).inner,par,'inner');
		filtering.(par.varnames{i}).inner.bciofull = logical(filtering.(par.varnames{i}).inner.filterMat'*(1*filtering.(par.varnames{i}).inner.bcio));
		
		[rhs.(par.varnames{i}).outer,filtering.(par.varnames{i}).outer.bcio] = bcf(grids,filtering,rhs.(par.varnames{i}).outer,par,'outer');
		filtering.(par.varnames{i}).outer.bciofull = logical(filtering.(par.varnames{i}).outer.filterMat'*(1*filtering.(par.varnames{i}).outer.bcio));
	
	end
	
end

function [rhs,bcio] = bcu(grids,filtering,rhs,par,side)
	
	flowrate = 1;
	
	on = filtering.u.(side).on;
	
	ymax = max(grids.u.(side).yinit);
	
	% driven flow
	drivy = ymax*ones(numel(on),1);
	in = ones(numel(on),flowrate);
	in(~(grids.u.(side).ymesh==drivy)) = 0;
	
	rhs = rhs + in;
	
	bcio = on&~on;
	
end

function [rhs,bcio] = bcv(grids,filtering,rhs,par,side)
	on = filtering.v.(side).on;
	bcio = on&~on;
end

function [rhs,bcio] = bcp(grids,filtering,rhs,par,side)
	on = filtering.p.(side).on;
	bcio = on&~on;
end

function [rhs,bcio] = bcq(grids,filtering,rhs,par,side)
	on = filtering.q.(side).on;
	bcio = on&~on;
end

