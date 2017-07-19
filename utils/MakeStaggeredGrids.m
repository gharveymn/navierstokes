function [grids,filtering,par] = MakeStaggeredGrids(par)
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
	
	%TODO just make a big rectangle and then cut it down
	
	%make dd bounds (if needed)
	if(par.ddrun)
		par.ddbounds{1}{1}(1) = xlimcoords(1);
		par.ddbounds{1}{2}(1) = xlimcoords(3) + par.ddoverlap;
		par.ddbounds{1}{1}(2) = ylimcoords(1);
		par.ddbounds{1}{2}(2) = ylimcoords(3);
		
		par.ddbounds{2}{1}(1) = xlimcoords(4);
		par.ddbounds{2}{2}(1) = par.h*round(1/par.h*(par.ddmidratio*xlimcoords(5)...
			+(1-par.ddmidratio)*xlimcoords(4))) + par.ddoverlap;
		par.ddbounds{2}{1}(2) = ylimcoords(7);
		par.ddbounds{2}{2}(2) = ylimcoords(4);
		
		par.ddbounds{3}{1}(1) = par.ddbounds{2}{2}(1) - par.ddoverlap;
		par.ddbounds{3}{2}(1) = xlimcoords(5);
		par.ddbounds{3}{1}(2) = ylimcoords(6);
		par.ddbounds{3}{2}(2) = ylimcoords(5);
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
	%hardcoded for symch
	for i=1:9
		if(i==5 || i==6)
			qxlimcoords(i) = qxlimcoords(i)-h;
		else
			qxlimcoords(i) = qxlimcoords(i)+h;
		end
		
		if(i>1 && i<6)
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
	%hardcoded for symch
	for i=1:9
		if(i==5 || i==6)
			pxlimcoords(i) = pxlimcoords(i)-h/2;
		else
			pxlimcoords(i) = pxlimcoords(i)+h/2;
		end
		
		if(i>1 && i<6)
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
	%hardcoded for symch
	for i=1:9
		if(i==5 || i==6)
			pxlimcoords(i) = pxlimcoords(i)+h/2;
		else
			pxlimcoords(i) = pxlimcoords(i)-h/2;
		end
		
		if(i>1 && i<6)
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
	%hardcoded for symch
	for i=1:9	
		if((i==5 || i==6))
			uxlimcoords(i) = uxlimcoords(i)-h;
		else
			uxlimcoords(i) = uxlimcoords(i)+h;
		end
		
		if(i>1 && i<6)
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
	for i=1:9
		if(i>1 && i<6)
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
	%hardcoded for symch
	for i=1:9
		if(i==5 || i==6)
			vxlimcoords(i) = vxlimcoords(i)-h/2;
		else
			vxlimcoords(i) = vxlimcoords(i)+h/2;
		end
		
		if(i>1 && i<6)
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
	%hardcoded for symch
	for i=1:9
		if(i==5 || i==6)
			vxlimcoords(i) = vxlimcoords(i)+h/2;
		else
			vxlimcoords(i) = vxlimcoords(i)-h/2;
		end
	end
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

function grids = placenums(grids,nx,ny)
	grids.nxp1 = nx+1;
	grids.nyp1 = ny+1;
	grids.nx = nx;
	grids.ny = ny;
	grids.nxm1 = nx-1;
	grids.nym1 = ny-1;
end
	
	
