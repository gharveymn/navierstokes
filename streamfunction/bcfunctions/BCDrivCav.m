function [rhs,filtering] = BCDrivCav(grids,filtering,rhs,par)
	
	[rhs.q.inner ,filtering] = bcq(grids,filtering,rhs.q.inner,par,'inner','inner');
	[rhs.q.outer1,filtering] = bcq(grids,filtering,rhs.q.outer1,par,'outer1','inner');
	[rhs.q.outer2,filtering] = bcq(grids,filtering,rhs.q.outer2,par,'outer2','inner');
	
end

function [rhs,filtering] = bcq(grids,filtering,rhs,par,side,sideinner)
	on = filtering.q.(side).on;
	
	xmin = min(grids.q.(side).xinit);
	xmax = max(grids.q.(side).xinit);	
	ymin = min(grids.q.(side).yinit);
	ymax = max(grids.q.(side).yinit);
	if(strcmp(side,'outer2'))
		%rhs(grids.q.outer2.ymesh==max(grids.q.outer2.ymesh)) = -1;
		rhs(grids.q.outer2.ymesh==ymax ...
			& ~(grids.q.outer2.xmesh==xmin | grids.q.outer2.xmesh==xmax)...
			& ~(grids.q.outer2.xmesh==min(grids.q.outer1.xinit) | grids.q.outer2.xmesh==max(grids.q.outer1.xinit)))...
			= -1;
	end
	filtering.q.(side).bcio = on&~on;
	filtering.q.(side).bciofull = logical(filtering.q.(side).filterMat'*(1*filtering.q.(side).bcio));
end

