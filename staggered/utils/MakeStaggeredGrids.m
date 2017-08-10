function [grids,filtering,par] = MakeStaggeredGrids(par)
	%MAKESTAGGEREDGRIDS parses the map file and gives us our mesh
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
	hx = par.hx;
	hy = par.hy;

	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	nx = (limits(2)-limits(1))/hx;
	ny = (limits(4)-limits(3))/hy;
	par.nx = nx;
	par.ny = ny;
	
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
	
	grids.u.inner.sz = struct('x',nx-1,'y',ny);
	grids.u.outer.sz = struct('x',nx+1,  'y',ny+2);
	
	grids.v.inner.sz = struct('x',nx,  'y',ny-1);
	grids.v.outer.sz = struct('x',nx+2,'y',ny+1);
	
	grids.p.inner.sz = struct('x',nx  ,'y',ny  );
	grids.p.outer.sz = struct('x',nx+2,'y',ny+2);
	
	grids.q.inner.sz = struct('x',nx-1,'y',ny-1);
	grids.q.outer.sz = struct('x',nx+1,'y',ny+1);

	filtering = struct;

	vns = par.varnames;
	for i=1:numel(vns)

		dx = (1 + par.nx - grids.(vns{i}).inner.sz.x)*hx/2;
		dy = (1 + par.ny - grids.(vns{i}).inner.sz.y)*hy/2;

		grids.(vns{i}).inner.xlimcoords = xlimcoords - dx*(incx - decx);
		grids.(vns{i}).inner.ylimcoords = ylimcoords - dy*(incy - decy);
		[grids,filtering] = CreateGrids(grids,filtering,par,'inner',vns{i});

		%outer
		dx = -(1 + par.nx - grids.(vns{i}).outer.sz.x)*hx/2;
		dy = -(1 + par.ny - grids.(vns{i}).outer.sz.y)*hy/2;

		grids.(vns{i}).outer.xlimcoords = xlimcoords + dx*(incx - decx);
		grids.(vns{i}).outer.ylimcoords = ylimcoords + dy*(incy - decy);
		[grids,filtering] = CreateGrids(grids,filtering,par,'outer',vns{i});

	end
	
	fclose('all');

end
	
	
