function [grids,filtering,par] = ParseValidIndices(par)
	%PARSEVALIDINDICES parses the map file and gives us our mesh
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
	
	%make streamfunction grid
	limits = [min(xlimcoords),max(xlimcoords),min(ylimcoords),max(ylimcoords)];
	qxinit = (limits(1):par.h:limits(2))';
	qyinit = (limits(3):par.h:limits(4))';
	nxp1 = numel(qxinit);
	nyp1 = numel(qyinit);
	[qgrids,qfiltering] = createGrids(qxinit,qyinit,nxp1,nyp1,xlimcoords,ylimcoords,h,true);
	
	%make pressure grids at cell centers
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
	
	[pgrids,pfiltering] = createGrids(pxinit,pyinit,nxp2,nyp2,pxlimcoords,pylimcoords,h,false);
	
	%make u velocity grid offset in the y direction
	uxinit = (limits(1):par.h:limits(2))';
	uyinit = (limits(3)-h/2:par.h:limits(4)+h/2)';
	nx = numel(uxinit);
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
	
	[ugrids,ufiltering] = createGrids(uxinit,uyinit,nx,nyp2,uxlimcoords,uylimcoords,h,false);
	
	%make u velocity grid offset in the y direction
	vxinit = (limits(1)-h/2:par.h:limits(2)+h/2)';
	vyinit = (limits(3):par.h:limits(4))';
	nxp2 = numel(vxinit);
	ny = numel(vyinit);
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
	[vgrids,vfiltering] = createGrids(vxinit,vyinit,nxp2,ny,vxlimcoords,vylimcoords,h,false);	
	
	
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
	
	grids = {qgrids,pgrids,ugrids,vgrids};
	filtering = {qfiltering,pfiltering,ufiltering,vfiltering};
	
	fclose('all');
end

function [grids,filtering] = createGrids(xinit,yinit,nx,ny,xlimcoords,ylimcoords,h,qflag)
	xmeshfull = kron(ones(ny,1),xinit);
	ymeshfull = kron(yinit,ones(nx,1));
	
	%Credit to Darren Engwirda for inpoly
	[valind,onfull] = inpoly(horzcat(xmeshfull,ymeshfull),horzcat(xlimcoords,ylimcoords));
	
	badcorners =		(effeq(xmeshfull,xlimcoords(1)) & effeq(ymeshfull,ylimcoords(1)))...
				|	(effeq(xmeshfull,xlimcoords(2)) & effeq(ymeshfull,ylimcoords(2)))...
				|	(effeq(xmeshfull,xlimcoords(4)) & effeq(ymeshfull,ylimcoords(4)))...
				|	(effeq(xmeshfull,xlimcoords(5)) & effeq(ymeshfull,ylimcoords(5)))...
				|	(effeq(xmeshfull,xlimcoords(6)) & effeq(ymeshfull,ylimcoords(6)))...
				|	(effeq(xmeshfull,xlimcoords(7)) & effeq(ymeshfull,ylimcoords(7)))...
				|	(effeq(xmeshfull,xlimcoords(9)) & effeq(ymeshfull,ylimcoords(9)));
	
	if(~qflag)
		valind = valind & ~badcorners;
		onfull = onfull & ~badcorners;
	end
	
	filterMat = spdiag(valind);
	filterMat = filterMat(valind,:);
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[nx,ny])';
	Ymesh = reshape(ymeshfull./valind,[nx,ny])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	grids = {xinit,yinit,xmesh,ymesh,Xmesh,Ymesh,xmeshfull,ymeshfull,nx,ny,h};
	
		
	%filterMat,{valindinner,valindouter},{on,onfull},{dbc,dbcfull},{gp1,gp2}
	filtering = {filterMat,{valind,valind},{on,onfull},{on,onfull},{onfull,[]}};
end

function bool=effeq(a,b)
	bool = abs(a-b) < eps;	
end
