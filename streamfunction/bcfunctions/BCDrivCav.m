function [rhs,filtering] = BCDrivCav(grids,filtering,rhs,par)
	
	[rhs.q.inner ,filtering] = bcq(grids,filtering,rhs.q.inner,par,'inner','inner');
	[rhs.q.outer1,filtering] = bcq(grids,filtering,rhs.q.outer1,par,'outer1','inner');
	[rhs.q.outer2,filtering] = bcq(grids,filtering,rhs.q.outer2,par,'outer2','inner');
	
end

function [rhs,filtering] = bcq(grids,filtering,rhs,par,side,sideinner)
	on = filtering.q.(side).on;
	if(strcmp(side,'outer2'))
		rhs(grids.q.outer2.ymesh==max(grids.q.outer2.ymesh)) = 1;
	end
	filtering.q.(side).bcio = on&~on;
	filtering.q.(side).bciofull = logical(filtering.q.(side).filterMat'*(1*filtering.q.(side).bcio));
end

