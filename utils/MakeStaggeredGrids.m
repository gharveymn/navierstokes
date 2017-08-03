function [grids,filtering,par] = MakeStaggeredGrids(par)
	%MakeGrids parses the map file and gives us our mesh
	%xinit,yinit are unmatched x and y sets -- vector
	%xmesh,ymesh have invalid indices removed -- vector
	%Xmesh,Ymesh have NaN wherever invalid -- matrix
	%filterMat filters out invalids when multiplied, inserts zeros when transposed and multiplied -- matrix
	%on holds the indices which are on the boundary -- in vector form, not prefiltered
	%
	%only availiable for SymCh
	
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
		[grids,filtering] = createGrids(grids,filtering,par,'inner',vns{i});

		%outer
		dx = -(1 + par.nx - grids.(vns{i}).outer.sz.x)*hx/2;
		dy = -(1 + par.ny - grids.(vns{i}).outer.sz.y)*hy/2;

		grids.(vns{i}).outer.xlimcoords = xlimcoords + dx*(incx - decx);
		grids.(vns{i}).outer.ylimcoords = ylimcoords + dy*(incy - decy);
		[grids,filtering] = createGrids(grids,filtering,par,'outer',vns{i});

	end
	
	fclose('all');

end

function [grids,filtering] = createGrids(grids,filtering,par,side,vn)
	
	hx = par.hx;
	hy = par.hy;
	
	szx = grids.(vn).(side).sz.x;
	szy = grids.(vn).(side).sz.y;

	xlimcoords = grids.(vn).(side).xlimcoords;
	ylimcoords = grids.(vn).(side).ylimcoords;

	
	xinit = linspace(min(grids.(vn).(side).xlimcoords),max(grids.(vn).(side).xlimcoords),szx)';
	yinit = linspace(min(grids.(vn).(side).ylimcoords),max(grids.(vn).(side).ylimcoords),szy)';
	
	xmeshfull = kron(ones(szy,1),xinit);
	ymeshfull = kron(yinit,ones(szx,1));
	
	[valind,onfull] = inpolygon(round(2*xmeshfull/hx),round(2*ymeshfull/hy),round(2*xlimcoords/hx),round(2*ylimcoords/hy));
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[szx,szy])';
	Ymesh = reshape(ymeshfull./valind,[szx,szy])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	grids.(vn).(side).xinit = xinit;
	grids.(vn).(side).yinit = yinit;
	grids.(vn).(side).xmesh = xmesh;
	grids.(vn).(side).ymesh = ymesh;
	grids.(vn).(side).Xmesh = Xmesh;
	grids.(vn).(side).Ymesh = Ymesh;
	grids.(vn).(side).xmeshfull = xmeshfull;
	grids.(vn).(side).ymeshfull = ymeshfull;
	
	filtering.(vn).(side).filterMat = filterMat;
	filtering.(vn).(side).valind = valind;
	filtering.(vn).(side).on = on;
	filtering.(vn).(side).onfull = onfull;
	
	[dbc,dbcfull] = boundarysides(grids.(vn),filtering.(vn),par,side,szx);
	
	filtering.(vn).(side).dbc = dbc;
	filtering.(vn).(side).dbcfull = dbcfull;
	filtering.(vn).(side).gp = [];

end
	
	
