function [grids,filtering,par] = MakeGrids(par)
	%MakeGrids parses the map file and gives us our mesh
	%xinit,yinit are unmatched x and y sets -- vector
	%xmesh,ymesh have invalid indices removed -- vector
	%Xmesh,Ymesh have NaN wherever invalid -- matrix
	%filterMat filters out invalids when multiplied, inserts zeros when transposed and multiplied -- matrix
	%on holds the indices which are on the boundary -- in vector form, not prefiltered
		
	file = fopen(par.mapfile, 'r');
	
	formatSpec = '%f';
	data = fscanf(file,formatSpec);
	
	%parse data into coordinates
	xlimcoords = zeros(numel(data)/2,1);
	ylimcoords = xlimcoords;
	
	for i=1:2:numel(data)
		xlimcoords((i+1)/2) = data(i);
		ylimcoords((i+1)/2) = data(i+1);
	end
	
	if(xlimcoords(1)~=xlimcoords(end) && ylimcoords(1)~=ylimcoords(end))
		xlimcoords = vertcat(xlimcoords,xlimcoords(1));
		ylimcoords = vertcat(ylimcoords,ylimcoords(1));
	end
	
	h = par.h;
	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	delx = limits(2)-limits(1);
	dely = limits(4)-limits(3);
	
	if(strcmp(par.griduse,'h'))
		hx = par.hx;
		hy = par.hy;
		nx = (delx)/hx;
		ny = (dely)/hy;
		par.nx = nx;
		par.ny = ny;
	else
		nx = par.nx;
		ny = par.ny;
		hx = (delx)/nx;
		hy = (dely)/ny;
		par.hx = hx;
		par.hy = hy;
	end
	
	c = convhull(xlimcoords,ylimcoords);
	convex = zeros(numel(xlimcoords),1);
	convex(c) = 1;
	convex(end) = convex(1);
	
	%find where we need to increase/decrease bounds
	[incx,onincx] = inpolygon(round((xlimcoords-hx)/hx),round(ylimcoords/hy),round(xlimcoords/hx),round(ylimcoords/hy));
	[decx,ondecx] = inpolygon(round((xlimcoords+hx)/hx),round(ylimcoords/hy),round(xlimcoords/hx),round(ylimcoords/hy));
	[incy,onincy] = inpolygon(round(xlimcoords/hx),round((ylimcoords-hy)/hy),round(xlimcoords/hx),round(ylimcoords/hy));
	[decy,ondecy] = inpolygon(round(xlimcoords/hx),round((ylimcoords+hy)/hy),round(xlimcoords/hx),round(ylimcoords/hy));
	
	incx = incx&(~onincx|convex);
	decx = decx&(~ondecx|convex);
	incy = incy&(~onincy|convex);
	decy = decy&(~ondecy|convex);
	
	grids.q.inner.sz = struct('x',nx-1,'y',ny-1);
	grids.q.outer1.sz = struct('x',nx+1,'y',ny+1);
	grids.q.outer2.sz = struct('x',nx+3,'y',ny+3);
	
	grids.q.hx = hx;
	grids.q.hy = hy;

	filtering = struct;
	
	dx = (1 + par.nx - grids.q.inner.sz.x)*hx/2;
	dy = (1 + par.ny - grids.q.inner.sz.y)*hy/2;

	grids.q.inner.xlimcoords = xlimcoords - dx*(incx - decx);
	grids.q.inner.ylimcoords = ylimcoords - dy*(incy - decy);
	[grids,filtering] = CreateGrids(grids,filtering,par,'inner','q');

	%outer
	dx = -(1 + par.nx - grids.q.outer1.sz.x)*hx/2;
	dy = -(1 + par.ny - grids.q.outer1.sz.y)*hy/2;

	grids.q.outer1.xlimcoords = xlimcoords + dx*(incx - decx);
	grids.q.outer1.ylimcoords = ylimcoords + dy*(incy - decy);
	[grids,filtering] = CreateGrids(grids,filtering,par,'outer1','q');
	
	dx = -(1 + par.nx - grids.q.outer2.sz.x)*hx/2;
	dy = -(1 + par.ny - grids.q.outer2.sz.y)*hy/2;
	
	grids.q.outer2.xlimcoords = xlimcoords + dx*(incx - decx);
	grids.q.outer2.ylimcoords = ylimcoords + dy*(incy - decy);
	[grids,filtering] = CreateGrids(grids,filtering,par,'outer2','q');
	
	
	fclose('all');
end