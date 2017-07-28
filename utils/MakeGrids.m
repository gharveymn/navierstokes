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
	[incx,onincx] = inpoly(horzcat(xlimcoords-h,ylimcoords),horzcat(xlimcoords,ylimcoords),[],h/4);
	[decx,ondecx] = inpoly(horzcat(xlimcoords+h,ylimcoords),horzcat(xlimcoords,ylimcoords),[],h/4);
	[incy,onincy] = inpoly(horzcat(xlimcoords,ylimcoords-h),horzcat(xlimcoords,ylimcoords),[],h/4);
	[decy,ondecy] = inpoly(horzcat(xlimcoords,ylimcoords+h),horzcat(xlimcoords,ylimcoords),[],h/4);
	
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
	[grids,filtering] = createGridsInner(qxinit,qyinit,nxm1,nym1,qxlimcoords,qylimcoords,h,par);

	%outer
	qxinit = (limits(1):par.h:limits(2))';
	qyinit = (limits(3):par.h:limits(4))';
	nxp1 = numel(qxinit);
	nyp1 = numel(qyinit);
	[grids,filtering] = createGridsOuter(qxinit,qyinit,nxp1,nyp1,xlimcoords,ylimcoords,h,grids,filtering,par);
	
	
	fclose('all');
end

function [grids,filtering] = createGridsInner(xinit,yinit,nx,ny,xlimcoords,ylimcoords,h,par)
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	%Credit to Darren Engwirda for inpoly
	[valind,onfull] = inpoly(horzcat(xmeshfull,ymeshfull),horzcat(xlimcoords,ylimcoords),[],h/4);
	
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
	[valind,onfull] = inpoly(horzcat(xmeshfull,ymeshfull),horzcat(xlimcoords,ylimcoords),[],h/4);
	
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