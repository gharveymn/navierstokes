function qmesh = NSFixY(grids,filtering,rhs,bc,par,solver)
	
	yinit = grids.yinit;
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	
	%for now we will just boost efficiency by keeping two types of derivative matrices

	nx = grids.nx;
	ny = grids.ny;
	h = grids.h;
	filterMat = filtering.filterMat;

	%make derivative matrices
	bih = biharmonic2(nx,ny,h,bc{1}{2}{1},bc{1}{2}{2});
	A = filterMat*bih*filterMat';
	
	%*debugging*
	%figure(1)
	%hold on
	
	for i=1:numel(yinit)
		ycurr = yinit(i);
		
		%*debugging*
		%scatter3(xmesh,ymesh,ymesh==ycurr,[],[ymesh==ycurr,zeros(numel(ymesh),1),1-(ymesh==ycurr)]);
		%drawnow;
		
		
		
		
	end
	
end
