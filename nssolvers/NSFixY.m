function psimesh = NSFixY(grids,filtering,rhs,bc,par,solver)
	
	yinit = grids{2};
	xmesh = grids{3};
	ymesh = grids{4};
	
	%for now we will just boost efficiency by keeping two types of derivative matrices

	nx = grids{9};
	ny = grids{10};
	h = grids{11};
	filterMat = filtering{1};

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

