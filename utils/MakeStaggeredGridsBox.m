function [grids,filtering,par] = MakeStaggeredGridsBox(par)
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
	
	
	h = par.h;

	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	par.nx = (limits(2)-limits(1))/h;
	par.ny = (limits(4)-limits(3))/h;
	par.nxp1 = par.nx+1;
	par.nyp1 = par.ny+1;
	
	
	%Q
	%------------------------------------------------
	
	%make streamfunction grid
	%inner
	qxinit = (limits(1)+h:par.h:limits(2)-h)';
	qyinit = (limits(3)+h:par.h:limits(4)-h)';
	nxm1 = numel(qxinit);
	nym1 = numel(qyinit);
	
	qxlimcoords = xlimcoords;
	qylimcoords = ylimcoords;
	%hardcoded for box
	for i=1:5
		if(i==2 || i==3)
			qxlimcoords(i) = qxlimcoords(i)-h;
		else
			qxlimcoords(i) = qxlimcoords(i)+h;
		end
		
		if(i>1 && i<4)
			qylimcoords(i) = qylimcoords(i)-h;
		else
			qylimcoords(i) = qylimcoords(i)+h;
		end
	end
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
	pxlimcoords = xlimcoords;
	pylimcoords = ylimcoords;
	%hardcoded for box
	for i=1:5
		if(i==2 || i==3)
			pxlimcoords(i) = pxlimcoords(i)-h/2;
		else
			pxlimcoords(i) = pxlimcoords(i)+h/2;
		end
		
		if(i>1 && i<4)
			pylimcoords(i) = pylimcoords(i)-h/2;
		else
			pylimcoords(i) = pylimcoords(i)+h/2;
		end
	end
	[pgrids,pfiltering] = createGridsInner(pxinit,pyinit,nx,ny,pxlimcoords,pylimcoords,h,par);

	%outer
	pxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	pyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nxp2 = numel(pxinit);
	nyp2 = numel(pyinit);
	pxlimcoords = xlimcoords;
	pylimcoords = ylimcoords;
	%hardcoded for box
	for i=1:5
		if(i==2 || i==3)
			pxlimcoords(i) = pxlimcoords(i)+h/2;
		else
			pxlimcoords(i) = pxlimcoords(i)-h/2;
		end
		
		if(i>1 && i<4)
			pylimcoords(i) = pylimcoords(i)+h/2;
		else
			pylimcoords(i) = pylimcoords(i)-h/2;
		end
	end
	[pgrids,pfiltering] = createGridsOuter(pxinit,pyinit,nxp2,nyp2,pxlimcoords,pylimcoords,h,pgrids,pfiltering,par);

	%U
	%------------------------------------------------

	%make u velocity grid offset in the y direction
	%inner
	uxinit = (limits(1)+h:par.h:limits(2)-h)';
	uyinit = (limits(3)+h/2:par.h:limits(4)-h/2)';
	nxm1 = numel(uxinit);
	ny = numel(uyinit);
	uxlimcoords = xlimcoords;
	uylimcoords = ylimcoords;
	%hardcoded for box
	for i=1:5	
		if((i==2 || i==3))
			uxlimcoords(i) = uxlimcoords(i)-h;
		else
			uxlimcoords(i) = uxlimcoords(i)+h;
		end
		
		if(i>1 && i<4)
			uylimcoords(i) = uylimcoords(i)-h/2;
		else
			uylimcoords(i) = uylimcoords(i)+h/2;
		end
	end
	[ugrids,ufiltering] = createGridsInner(uxinit,uyinit,nxm1,ny,uxlimcoords,uylimcoords,h,par);

	%outer
	uxinit = (limits(1):par.h:limits(2))';
	uyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nxp1 = numel(uxinit);
	nyp2 = numel(uyinit);
	uxlimcoords = xlimcoords;
	uylimcoords = ylimcoords;
	%hardcoded for symch
	for i=1:5
		if(i>1 && i<4)
			uylimcoords(i) = uylimcoords(i)+h/2;
		else
			uylimcoords(i) = uylimcoords(i)-h/2;
		end
	end
	
	[ugrids,ufiltering] = createGridsOuter(uxinit,uyinit,nxp1,nyp2,uxlimcoords,uylimcoords,h,ugrids,ufiltering,par);
	
	%V
	%------------------------------------------------

	%make v velocity grid offset in the x direction
	%inner
	vxinit = (limits(1)+h/2:par.h:limits(2)-h/2)';
	vyinit = (limits(3)+h:par.h:limits(4)-h)';
	nx = numel(vxinit);
	nym1 = numel(vyinit);
	vxlimcoords = xlimcoords;
	vylimcoords = ylimcoords;
	%hardcoded for box
	for i=1:5
		if(i==2 || i==3)
			vxlimcoords(i) = vxlimcoords(i)-h/2;
		else
			vxlimcoords(i) = vxlimcoords(i)+h/2;
		end
		
		if(i>1 && i<4)
			vylimcoords(i) = vylimcoords(i)-h;
		else
			vylimcoords(i) = vylimcoords(i)+h;
		end
	end
	[vgrids,vfiltering] = createGridsInner(vxinit,vyinit,nx,nym1,vxlimcoords,vylimcoords,h,par);	

	%outer
	vxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	vyinit = (limits(3):par.h:limits(4))';
	nxp2 = numel(vxinit);
	nyp1 = numel(vyinit);
	vxlimcoords = xlimcoords;
	vylimcoords = ylimcoords;
	%hardcoded for box
	for i=1:5
		if(i==2 || i==3)
			vxlimcoords(i) = vxlimcoords(i)+h/2;
		else
			vxlimcoords(i) = vxlimcoords(i)-h/2;
		end
	end
	[vgrids,vfiltering] = createGridsOuter(vxinit,vyinit,nxp2,nyp1,vxlimcoords,vylimcoords,h,vgrids,vfiltering,par);	

	%------------------------------------------------
	
	
	%NOTE: filtering{4} ie {bc,bcfull} will change to become 2x2 dimensional
	%		this is a hack because I really do not want to touch closure since that is also a hack
	
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

	xmin = min(xinit);
	xmax = max(xinit);
	ymin = min(yinit);
	ymax = max(yinit);
	
	valind = logical(ones(numel(xmeshfull),1));
	onfull = xmeshfull==xmin|xmeshfull==xmax|ymeshfull==ymin|ymeshfull==ymax;
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull;
	
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
	
	dbc.w = xmesh==xmin;
	dbc.e = xmesh==xmax;
	dbc.s = ymesh==ymin;
	dbc.n = ymesh==ymax;
	dbc.c = (dbc.w|dbc.e)&(dbc.s|dbc.n);

	dbcfull.w = xmeshfull==xmin;
	dbcfull.e = xmeshfull==xmax;
	dbcfull.s = ymeshfull==ymin;
	dbcfull.n = ymeshfull==ymax;
	dbcfull.c = (dbcfull.w|dbcfull.e)&(dbcfull.s|dbcfull.n);
	
	filtering.inner.dbc = dbc;
	filtering.inner.dbcfull = dbcfull;
	filtering.inner.gp = [];

end

function [grids,filtering] = createGridsOuter(xinit,yinit,nx,ny,xlimcoords,ylimcoords,h,grids,filtering,par)
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	xmin = min(xinit);
	xmax = max(xinit);
	ymin = min(yinit);
	ymax = max(yinit);
	
	valind = logical(ones(numel(xmeshfull),1));
	onfull = xmeshfull==xmin|xmeshfull==xmax|ymeshfull==ymin|ymeshfull==ymax;
	
	% if(~qflag)
	% 	valind = valind & ~badcorners;
	% 	onfull = onfull & ~badcorners;
	% end
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull;
	
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
	
	dbc.w = xmesh==xmin;
	dbc.e = xmesh==xmax;
	dbc.s = ymesh==ymin;
	dbc.n = ymesh==ymax;
	dbc.c = (dbc.w|dbc.e)&(dbc.s|dbc.n);

	dbcfull.w = xmeshfull==xmin;
	dbcfull.e = xmeshfull==xmax;
	dbcfull.s = ymeshfull==ymin;
	dbcfull.n = ymeshfull==ymax;
	dbcfull.c = (dbcfull.w|dbcfull.e)&(dbcfull.s|dbcfull.n);

	
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
	
	
