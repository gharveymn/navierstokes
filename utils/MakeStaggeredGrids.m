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
	
	%Q
	%------------------------------------------------
	
	%make streamfunction grid
	%inner
	qxlimcoords = xlimcoords - hx*incx + hx*decx;
	qylimcoords = ylimcoords - hy*incy + hy*decy;
	[qgrids,qfiltering] = createGrids(limits(1)+hx,limits(2)-hx,limits(3)+hy,limits(4)-hy,nx-1,ny-1,qxlimcoords,qylimcoords,hx,hy,par,'inner');	

	%outer
	[qgrids,qfiltering] = createGrids(limits(1),limits(2),limits(3),limits(4),nx+1,ny+1,xlimcoords,ylimcoords,hx,hy,par,'outer',qgrids,qfiltering);	
	
	%P
	%------------------------------------------------

	%make pressure grids at cell centers
	%inner
	pxlimcoords = xlimcoords - hx/2*incx + hx/2*decx;
	pylimcoords = ylimcoords - hy/2*incy + hy/2*decy;
	[pgrids,pfiltering] = createGrids(limits(1)+hx/2,limits(2)-hx/2,limits(3)+hy/2,limits(4)-hy/2,nx,ny,pxlimcoords,pylimcoords,hx,hy,par,'inner');

	%outer
	pxlimcoords = xlimcoords + hx/2*incx - hx/2*decx;
	pylimcoords = ylimcoords + hy/2*incy - hy/2*decy;
	[pgrids,pfiltering] = createGrids(limits(1)-hx/2,limits(2)+hx/2,limits(3)-hy/2,limits(4)+hy/2,nx+2,ny+2,pxlimcoords,pylimcoords,hx,hy,par,'outer',pgrids,pfiltering);	

	%U
	%------------------------------------------------

	%make u velocity grid offset in the y direction
	%inner
	uxlimcoords = xlimcoords - hx*incx + hx*decx;
	uylimcoords = ylimcoords - hy/2*incy + hy/2*decy;
	[ugrids,ufiltering] = createGrids(limits(1)+hx,limits(2)-hx,limits(3)+hy/2,limits(4)-hy/2,nx-1,ny,uxlimcoords,uylimcoords,hx,hy,par,'inner');

	%outer
	uxlimcoords = xlimcoords;
	uylimcoords = ylimcoords + hy/2*incy - hy/2*decy;	
	[ugrids,ufiltering] = createGrids(limits(1),limits(2),limits(3)-hy/2,limits(4)+hy/2,nx+1,ny+2,uxlimcoords,uylimcoords,hx,hy,par,'outer',ugrids,ufiltering);	
	
	%V
	%------------------------------------------------

	%make v velocity grid offset in the x direction
	%inner
	vxlimcoords = xlimcoords - hx/2*incx + hx/2*decx;
	vylimcoords = ylimcoords - hy*incy + hy*decy;
	[vgrids,vfiltering] = createGrids(limits(1)+hx/2,limits(2)-hx/2,limits(3)+hy,limits(4)-hy,nx,ny-1,vxlimcoords,vylimcoords,hx,hy,par,'inner');	

	%outer
	vxlimcoords = xlimcoords + hx/2*incx - hx/2*decx;
	vylimcoords = ylimcoords;
	[vgrids,vfiltering] = createGrids(limits(1)-hx/2,limits(2)+hx/2,limits(3),limits(4),nx+2,ny+1,vxlimcoords,vylimcoords,hx,hy,par,'outer',vgrids,vfiltering);
	
	grids.p = pgrids;
	grids.u = ugrids;
	grids.v = vgrids;
	grids.q = qgrids;
	filtering.p = pfiltering;
	filtering.u = ufiltering;
	filtering.v = vfiltering;
	filtering.q = qfiltering;
	
	%this is absolutely spaghetti code at this point
	grids.p = placenums(grids.p,nx,ny);
	grids.u = placenums(grids.u,nx,ny);
	grids.v = placenums(grids.v,nx,ny);
	grids.q = placenums(grids.q,nx,ny);
	
	fclose('all');
end

function [grids,filtering] = createGrids(x1,x2,y1,y2,nx,ny,xlimcoords,ylimcoords,hx,hy,par,side,grids,filtering)
	
	xinit = linspace(x1,x2,nx)';
	yinit = linspace(y1,y2,ny)';
	
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	a = round(2*xmeshfull/hx);
	b = round(2*ymeshfull/hy);
	c = round(2*xlimcoords/hx);
	d = round(2*ylimcoords/hy);
	
	[valind,onfull] = inpolygon(round(2*xmeshfull/hx),round(2*ymeshfull/hy),round(2*xlimcoords/hx),round(2*ylimcoords/hy));
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[nx,ny])';
	Ymesh = reshape(ymeshfull./valind,[nx,ny])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	g.xinit = xinit;
	g.yinit = yinit;
	g.xmesh = xmesh;
	g.ymesh = ymesh;
	g.Xmesh = Xmesh;
	g.Ymesh = Ymesh;
	g.xmeshfull = xmeshfull;
	g.ymeshfull = ymeshfull;
	grids.(side) = g;
	
	f.filterMat = filterMat;
	f.valind = valind;
	f.on = on;
	f.onfull = onfull;
	filtering.(side) = f;
	
	[dbc,dbcfull] = boundarysides(grids,filtering,par,side,nx);
	
	filtering.(side).dbc = dbc;
	filtering.(side).dbcfull = dbcfull;
	filtering.(side).gp = [];

end

function grids = placenums(grids,nx,ny)
	grids.nxp1 = nx+1;
	grids.nyp1 = ny+1;
	grids.nx = nx;
	grids.ny = ny;
	grids.nxm1 = nx-1;
	grids.nym1 = ny-1;
end
	
	
