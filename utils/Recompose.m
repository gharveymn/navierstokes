function qmesh = Recompose(grids,newqs,par)
	%RECOMPOSE for use with symch
	
	ddbounds = par.ddbounds;
	del = (par.order-1)*par.h+par.h/2;
	
	%just to allocate
	qmesh = 0*grids.xmesh;
	
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	
	b1 = ddbounds{1};
	b2 = ddbounds{2};
	b3 = ddbounds{3};
	
	%maybe we'll average out overlap, but convergence should take care of that
	qmesh((xmesh >= b1{1}(1)-del) & (xmesh <= b1{2}(1)+del) & (ymesh >= b1{1}(2)-del) & (ymesh <= b1{2}(2)+del)) = newqs{1};
	qmesh((xmesh >= b2{1}(1)-del) & (xmesh <= b2{2}(1)+del) & (ymesh >= b2{1}(2)-del) & (ymesh <= b2{2}(2)+del)) = newqs{2};
	qmesh((xmesh >= b3{1}(1)-del) & (xmesh <= b3{2}(1)+del) & (ymesh >= b3{1}(2)-del) & (ymesh <= b3{2}(2)+del)) = newqs{3};

end

