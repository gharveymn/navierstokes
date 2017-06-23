function psimesh = Recompose(grids,newpsis,par)
	%RECOMPOSE for use with symch
	
	ddbounds = par.ddbounds;
	del = (par.order-1)*par.h+par.h/2;
	
	%just to allocate
	psimesh = 0*grids{3};
	
	xmesh = grids{3};
	ymesh = grids{4};
	
	b1 = ddbounds{1};
	b2 = ddbounds{2};
	b3 = ddbounds{3};
	
	%maybe we'll average out overlap, but convergence should take care of that
	psimesh((xmesh >= b1{1}(1)-del) & (xmesh <= b1{2}(1)+del) & (ymesh >= b1{1}(2)-del) & (ymesh <= b1{2}(2)+del)) = newpsis{1};
	psimesh((xmesh >= b2{1}(1)-del) & (xmesh <= b2{2}(1)+del) & (ymesh >= b2{1}(2)-del) & (ymesh <= b2{2}(2)+del)) = newpsis{2};
	psimesh((xmesh >= b3{1}(1)-del) & (xmesh <= b3{2}(1)+del) & (ymesh >= b3{1}(2)-del) & (ymesh <= b3{2}(2)+del)) = newpsis{3};

end

