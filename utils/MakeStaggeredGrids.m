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

	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	par.nx = (limits(2)-limits(1))/h;
	par.ny = (limits(4)-limits(3))/h;
	par.nxp1 = par.nx+1;
	par.nyp1 = par.ny+1;
	
	c = convhull(xlimcoords,ylimcoords);
	convex = zeros(numel(xlimcoords),1);
	convex(c) = 1;
	convex(end) = convex(1);
	
	%find where we need to increase/decrease bounds
	[incx,onincx] = inpolygon(round((xlimcoords-h)/h),round(ylimcoords/h),round(xlimcoords/h),round(ylimcoords/h));
	[decx,ondecx] = inpolygon(round((xlimcoords+h)/h),round(ylimcoords/h),round(xlimcoords/h),round(ylimcoords/h));
	[incy,onincy] = inpolygon(round(xlimcoords/h),round((ylimcoords-h)/h),round(xlimcoords/h),round(ylimcoords/h));
	[decy,ondecy] = inpolygon(round(xlimcoords/h),round((ylimcoords+h)/h),round(xlimcoords/h),round(ylimcoords/h));
	
	incx = incx&(~onincx|convex);
	decx = decx&(~ondecx|convex);
	incy = incy&(~onincy|convex);
	decy = decy&(~ondecy|convex);
	
	%Q
	%------------------------------------------------
	
	%make streamfunction grid
	%inner
	qxinit = (limits(1)+h:par.h:limits(2)-h)';
	qyinit = (limits(3)+h:par.h:limits(4)-h)';
	nxm1 = numel(qxinit);
	nym1 = numel(qyinit);
	
	qxlimcoords = xlimcoords - h*incx + h*decx;
	qylimcoords = ylimcoords - h*incy + h*decy;
	[qgrids,qfiltering] = createGrids(qxinit,qyinit,nxm1,nym1,qxlimcoords,qylimcoords,h,par,'inner');	

	%outer
	qxinit = (limits(1):par.h:limits(2))';
	qyinit = (limits(3):par.h:limits(4))';
	nxp1 = numel(qxinit);
	nyp1 = numel(qyinit);
	[qgrids,qfiltering] = createGrids(qxinit,qyinit,nxp1,nyp1,xlimcoords,ylimcoords,h,par,'outer',qgrids,qfiltering);	
	
	%P
	%------------------------------------------------

	%make pressure grids at cell centers
	%inner
	pxinit = (limits(1)+h/2:par.h:limits(2)-h/2)';
	pyinit = (limits(3)+h/2:par.h:limits(4)-h/2)';
	nx = numel(pxinit);
	ny = numel(pyinit);
	pxlimcoords = xlimcoords - h/2*incx + h/2*decx;
	pylimcoords = ylimcoords - h/2*incy + h/2*decy;
	[pgrids,pfiltering] = createGrids(pxinit,pyinit,nx,ny,pxlimcoords,pylimcoords,h,par,'inner');

	%outer
	pxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	pyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nxp2 = numel(pxinit);
	nyp2 = numel(pyinit);
	pxlimcoords = xlimcoords + h/2*incx - h/2*decx;
	pylimcoords = ylimcoords + h/2*incy - h/2*decy;
	[pgrids,pfiltering] = createGrids(pxinit,pyinit,nxp2,nyp2,pxlimcoords,pylimcoords,h,par,'outer',pgrids,pfiltering);	

	%U
	%------------------------------------------------

	%make u velocity grid offset in the y direction
	%inner
	uxinit = (limits(1)+h:par.h:limits(2)-h)';
	uyinit = (limits(3)+h/2:par.h:limits(4)-h/2)';
	nxm1 = numel(uxinit);
	ny = numel(uyinit);
	uxlimcoords = xlimcoords - h*incx + h*decx;
	uylimcoords = ylimcoords - h/2*incy + h/2*decy;
	[ugrids,ufiltering] = createGrids(uxinit,uyinit,nxm1,ny,uxlimcoords,uylimcoords,h,par,'inner');

	%outer
	uxinit = (limits(1):par.h:limits(2))';
	uyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nxp1 = numel(uxinit);
	nyp2 = numel(uyinit);
	uxlimcoords = xlimcoords;
	uylimcoords = ylimcoords + h/2*incy - h/2*decy;	
	[ugrids,ufiltering] = createGrids(uxinit,uyinit,nxp1,nyp2,uxlimcoords,uylimcoords,h,par,'outer',ugrids,ufiltering);	
	
	%V
	%------------------------------------------------

	%make v velocity grid offset in the x direction
	%inner
	vxinit = (limits(1)+h/2:par.h:limits(2)-h/2)';
	vyinit = (limits(3)+h:par.h:limits(4)-h)';
	nx = numel(vxinit);
	nym1 = numel(vyinit);
	vxlimcoords = xlimcoords - h/2*incx + h/2*decx;
	vylimcoords = ylimcoords - h*incy + h*decy;
	[vgrids,vfiltering] = createGrids(vxinit,vyinit,nx,nym1,vxlimcoords,vylimcoords,h,par,'inner');	

	%outer
	vxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	vyinit = (limits(3):par.h:limits(4))';
	nxp2 = numel(vxinit);
	nyp1 = numel(vyinit);
	vxlimcoords = xlimcoords + h/2*incx - h/2*decx;
	vylimcoords = ylimcoords;
	[vgrids,vfiltering] = createGrids(vxinit,vyinit,nxp2,nyp1,vxlimcoords,vylimcoords,h,par,'outer',vgrids,vfiltering);
	
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

function [grids,filtering] = createGrids(xinit,yinit,nx,ny,xlimcoords,ylimcoords,h,par,side,grids,filtering)
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	[valind,onfull] = inpolygon(round(2*xmeshfull/h),round(2*ymeshfull/h),round(2*xlimcoords/h),round(2*ylimcoords/h));
	
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
	grids.h = h;
	
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
	
	
