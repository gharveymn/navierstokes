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
	
	vns = par.varnames;
	
	h = par.h;
	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	delx = limits(2)-limits(1);
	dely = limits(4)-limits(3);
	
	if(strcmp(par.griduse,'h'))
		hx = par.hx;
		hy = par.hy;
		nx = delx/hx;
		ny = dely/hy;
		par.nx = nx;
		par.ny = ny;
	else
		nx = par.nx;
		ny = par.ny;
		hx = delx/nx;
		hy = dely/ny;
		par.hx = hx;
		par.hy = hy;
	end
	
	if(par.useinterp)
		pnx = par.pnx;
		pny = par.pny;
		phx = delx/pnx;
		phy = dely/pny;
		par.phx = phx;
		par.phy = phy;		
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
	
	grids.u.inner.sz = struct('x',nx-1,'y',ny);
	grids.u.outer.sz = struct('x',nx+1,  'y',ny+2);
	
	grids.v.inner.sz = struct('x',nx,  'y',ny-1);
	grids.v.outer.sz = struct('x',nx+2,'y',ny+1);
	
	grids.p.inner.sz = struct('x',nx  ,'y',ny  );
	grids.p.outer.sz = struct('x',nx+2,'y',ny+2);
	
	grids.q.inner.sz = struct('x',nx-1,'y',ny-1);
	grids.q.outer.sz = struct('x',nx+1,'y',ny+1);
	
	for i=1:numel(vns)
		grids.(vns{i}).hx = hx;
		grids.(vns{i}).hy = hy;
	end

	filtering = struct;

	for i=1:numel(vns)
		vn = vns{i};

		dx = (1 + nx - grids.(vn).inner.sz.x)*hx/2;
		dy = (1 + ny - grids.(vn).inner.sz.y)*hy/2;

		grids.(vn).inner.xlimcoords = xlimcoords - dx*(incx - decx);
		grids.(vn).inner.ylimcoords = ylimcoords - dy*(incy - decy);
		[grids,filtering] = CreateGrids(grids,filtering,par,'inner',vns{i});

		%outer
		dx = -(1 + nx - grids.(vns{i}).outer.sz.x)*hx/2;
		dy = -(1 + ny - grids.(vns{i}).outer.sz.y)*hy/2;

		grids.(vn).outer.xlimcoords = xlimcoords + dx*(incx - decx);
		grids.(vn).outer.ylimcoords = ylimcoords + dy*(incy - decy);
		[grids,filtering] = CreateGrids(grids,filtering,par,'outer',vns{i});

	end
	
	if(par.useinterp)
		
		gplot.u.inner.sz = struct('x',pnx-1,'y',pny);
		gplot.u.outer.sz = struct('x',pnx+1,  'y',pny+2);

		gplot.v.inner.sz = struct('x',pnx,  'y',pny-1);
		gplot.v.outer.sz = struct('x',pnx+2,'y',pny+1);

		gplot.p.inner.sz = struct('x',pnx  ,'y',pny  );
		gplot.p.outer.sz = struct('x',pnx+2,'y',pny+2);

		gplot.q.inner.sz = struct('x',pnx-1,'y',pny-1);
		gplot.q.outer.sz = struct('x',pnx+1,'y',pny+1);
		
		for i=1:numel(vns)
			gplot.(vns{i}).hx = phx;
			gplot.(vns{i}).hy = phy;
		end
		
		fplot = struct;
		
		for i=1:numel(vns)
			
			vn = vns{i};

			dx = (1 + pnx - gplot.(vn).inner.sz.x)*phx/2;
			dy = (1 + pny - gplot.(vn).inner.sz.y)*phy/2;

			gplot.(vn).inner.xlimcoords = xlimcoords - dx*(incx - decx);
			gplot.(vn).inner.ylimcoords = ylimcoords - dy*(incy - decy);
			[gplot,fplot] = CreateGrids(gplot,fplot,par,'inner',vn);

			%outer
			dx = -(1 + pnx - gplot.(vn).outer.sz.x)*phx/2;
			dy = -(1 + pny - gplot.(vn).outer.sz.y)*phy/2;

			gplot.(vn).outer.xlimcoords = xlimcoords + dx*(incx - decx);
			gplot.(vn).outer.ylimcoords = ylimcoords + dy*(incy - decy);
			[gplot,fplot] = CreateGrids(gplot,fplot,par,'outer',vn);

		end
		
		grids.plot = gplot;
		filtering.plot = fplot;
		
	end
	
	fclose('all');

end
	
	
