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
	[incx,onincx] = inpolygon(xlimcoords-h,ylimcoords,xlimcoords,ylimcoords);
	[decx,ondecx] = inpolygon(xlimcoords+h,ylimcoords,xlimcoords,ylimcoords);
	[incy,onincy] = inpolygon(xlimcoords,ylimcoords-h,xlimcoords,ylimcoords);
	[decy,ondecy] = inpolygon(xlimcoords,ylimcoords+h,xlimcoords,ylimcoords);
	
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
	[qgrids,qfiltering] = createGridsInner(qxinit,qyinit,nxm1,nym1,qxlimcoords,qylimcoords,h,par);

	%outer
	qxinit = (limits(1):par.h:limits(2))';
	qyinit = (limits(3):par.h:limits(4))';
	nxp1 = numel(qxinit);
	nyp1 = numel(qyinit);
	[qgrids,qfiltering] = createGridsOuter(qxinit,qyinit,nxp1,nyp1,xlimcoords,ylimcoords,h,qgrids,qfiltering,par);
	
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
	[pgrids,pfiltering] = createGridsInner(pxinit,pyinit,nx,ny,pxlimcoords,pylimcoords,h,par);

	%outer
	pxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	pyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nxp2 = numel(pxinit);
	nyp2 = numel(pyinit);
	pxlimcoords = xlimcoords + h/2*incx - h/2*decx;
	pylimcoords = ylimcoords + h/2*incy - h/2*decy;
	[pgrids,pfiltering] = createGridsOuter(pxinit,pyinit,nxp2,nyp2,pxlimcoords,pylimcoords,h,pgrids,pfiltering,par);

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
	[ugrids,ufiltering] = createGridsInner(uxinit,uyinit,nxm1,ny,uxlimcoords,uylimcoords,h,par);

	%outer
	uxinit = (limits(1):par.h:limits(2))';
	uyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nxp1 = numel(uxinit);
	nyp2 = numel(uyinit);
	uxlimcoords = xlimcoords;
	uylimcoords = ylimcoords + h/2*incy - h/2*decy;	
	[ugrids,ufiltering] = createGridsOuter(uxinit,uyinit,nxp1,nyp2,uxlimcoords,uylimcoords,h,ugrids,ufiltering,par);
	
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
	[vgrids,vfiltering] = createGridsInner(vxinit,vyinit,nx,nym1,vxlimcoords,vylimcoords,h,par);	

	%outer
	vxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	vyinit = (limits(3):par.h:limits(4))';
	nxp2 = numel(vxinit);
	nyp1 = numel(vyinit);
	vxlimcoords = xlimcoords + h/2*incx - h/2*decx;
	vylimcoords = ylimcoords;
	[vgrids,vfiltering] = createGridsOuter(vxinit,vyinit,nxp2,nyp1,vxlimcoords,vylimcoords,h,vgrids,vfiltering,par);	

	%------------------------------------------------
	
	
	%NOTE: filtering{4} ie {bc,bcfull} will change to become 2x2 dimensional
	%		this is a hack because I really do not want to touch closure since that is also a hack
	
	%we did some arithmetic up there so just
	%make sure ddbounds are actually in the grids
	[~,in] = min(abs(qxinit-par.ddbounds{1}{2}(1)));
	par.ddbounds{1}{2}(1) = qxinit(in);
	
	[~,in] = min(abs(qxinit-par.ddbounds{2}{2}(1)));
	par.ddbounds{2}{2}(1) = qxinit(in);
	
	[~,in] = min(abs(qxinit-par.ddbounds{3}{1}(1)));
	par.ddbounds{3}{1}(1) = qxinit(in);
	
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

function [grids,filtering] = createGridsInner(xinit,yinit,nx,ny,xlimcoords,ylimcoords,h,par)
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	%Credit to Darren Engwirda for inpoly
	[valind,onfull] = inpolygon(xmeshfull,ymeshfull,xlimcoords,ylimcoords);
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[nx,ny])';
	Ymesh = reshape(ymeshfull./valind,[nx,ny])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	inner.xinit = xinit;
	inner.yinit = yinit;
	inner.xmesh = xmesh;
	inner.ymesh = ymesh;
	inner.Xmesh = Xmesh;
	inner.Ymesh = Ymesh;
	inner.xmeshfull = xmeshfull;
	inner.ymeshfull = ymeshfull;
	grids.inner = inner;
	grids.h = h;
	
	finner.filterMat = filterMat;
	finner.valind = valind;
	finner.on = on;
	finner.onfull = onfull;
	filtering.inner = finner;
	
	[dbc,dbcfull] = boundarysides(grids,filtering,par,'inner',nx);
	
	filtering.inner.dbc = dbc;
	filtering.inner.dbcfull = dbcfull;
	filtering.inner.gp = [];

end

function [grids,filtering] = createGridsOuter(xinit,yinit,nx,ny,xlimcoords,ylimcoords,h,grids,filtering,par)
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	%Credit to Darren Engwirda for inpoly
	[valind,onfull] = inpolygon(xmeshfull,ymeshfull,xlimcoords,ylimcoords);
	
	% badcorners =		(effeq(xmeshfull,xlimcoords(1)) & effeq(ymeshfull,ylimcoords(1)))...
	% 			|	(effeq(xmeshfull,xlimcoords(2)) & effeq(ymeshfull,ylimcoords(2)))...
	% 			|	(effeq(xmeshfull,xlimcoords(4)) & effeq(ymeshfull,ylimcoords(4)))...
	% 			|	(effeq(xmeshfull,xlimcoords(5)) & effeq(ymeshfull,ylimcoords(5)))...
	% 			|	(effeq(xmeshfull,xlimcoords(6)) & effeq(ymeshfull,ylimcoords(6)))...
	% 			|	(effeq(xmeshfull,xlimcoords(7)) & effeq(ymeshfull,ylimcoords(7)))...
	% 			|	(effeq(xmeshfull,xlimcoords(9)) & effeq(ymeshfull,ylimcoords(9)));
	
	% if(~qflag)
	% 	valind = valind & ~badcorners;
	% 	onfull = onfull & ~badcorners;
	% end
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[nx,ny])';
	Ymesh = reshape(ymeshfull./valind,[nx,ny])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	outer.xinit = xinit;
	outer.yinit = yinit;
	outer.xmesh = xmesh;
	outer.ymesh = ymesh;
	outer.Xmesh = Xmesh;
	outer.Ymesh = Ymesh;
	outer.xmeshfull = xmeshfull;
	outer.ymeshfull = ymeshfull;
	grids.outer = outer;

	grids.h = h;
	
	fouter.filterMat = filterMat;
	fouter.valind = valind;
	fouter.on = on;
	fouter.onfull = onfull;
	filtering.outer = fouter;
	
	[dbc,dbcfull] = boundarysides(grids,filtering,par,'outer',nx);
	
	filtering.outer.dbc = dbc;
	filtering.outer.dbcfull = dbcfull;
	filtering.outer.gp = [];

end

function grids = placenums(grids,nx,ny)
	grids.nxp1 = nx+1;
	grids.nyp1 = ny+1;
	grids.nx = nx;
	grids.ny = ny;
	grids.nxm1 = nx-1;
	grids.nym1 = ny-1;
end
	
	
