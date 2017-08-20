function [grids,filtering] = CreateGrids(grids,filtering,par,side,vn)
	
	hx = grids.(vn).hx;
	hy = grids.(vn).hy;
	
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
	
	valmat = reshape(filtering.u.inner.valind,[szx,szy])';
	
	on = onfull(valind);
	
	Xmesh = reshape(xmeshfull./valind,[szx,szy])';
	Ymesh = reshape(ymeshfull./valind,[szx,szy])';
	
	xmesh = filterMat*xmeshfull;
	ymesh = filterMat*ymeshfull;
	
	if(par.useinterp)
		
	end
	
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
	filtering.(vn).(side).valmat = valmat;
	filtering.(vn).(side).on = on;
	filtering.(vn).(side).onfull = onfull;
	
	[dbc,dbcfull] = boundarysides(grids.(vn),filtering.(vn),par,side,szx);
	
	filtering.(vn).(side).dbc = dbc;
	filtering.(vn).(side).dbcfull = dbcfull;
	filtering.(vn).(side).gp = [];

end