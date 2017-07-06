function [dbc,dbcfull] = boundarysides(grids,filtering,par,side)
	%GETWHEREBOUNDARIES I'm somewhat suprised this actually works
	
	if(~exist('side','var'))
		side = 'inner';
	end
	
	if(strcmp(side,'inner'))
		xmeshfull = grids.inner.xmeshfull;
		ymeshfull = grids.inner.ymeshfull;
		onfull = filtering.inner.onfull;
		valindinner = filtering.inner.valind;
	elseif(strcmp(side,'outer'))
		xmeshfull = grids.outer.xmeshfull;
		ymeshfull = grids.outer.ymeshfull;
		onfull = filtering.outer.onfull;
		valindouter = filtering.outer.valind;
	else
		ME = MException('boundarysides:invalidParameterException','Invalid value for side');
		throw(ME)
	end
	nxp1 = par.nxp1;
	
	xmin = min(xmeshfull);
	xmax = max(xmeshfull);
	ymin = min(ymeshfull);
	ymax = max(ymeshfull);
	
	%set to contain all possible indices
	bcw = onfull;
	bce = onfull;
	bcs = onfull;
	bcn = onfull;
	
	%get indices on the ends
	xminb = (xmeshfull==xmin);
	xmaxb = (xmeshfull==xmax);
	yminb = (ymeshfull==ymin);
	ymaxb = (ymeshfull==ymax);
	
	%if an index is on the overall boundary then it must be a boundary point
	bcw = bcw&xminb;
	bce = bce&xmaxb;
	bcs = bcs&yminb;
	bcn = bcn&ymaxb;
	
	%make shifted index matrices, filter out those indices at the max values since circshift loops
	
	if(strcmp(side,'inner'))
		r = circshift(valindinner&~xmaxb,1);
		l = circshift(valindinner&~xminb,-1);
		u = circshift(valindinner&~ymaxb,nxp1);
		d = circshift(valindinner&~yminb,-nxp1);
	else
		r = circshift(valindouter&~xmaxb,1);
		l = circshift(valindouter&~xminb,-1);
		u = circshift(valindouter&~ymaxb,nxp1);
		d = circshift(valindouter&~yminb,-nxp1);
	end
	
	%if the shift makes it go off the boundary then we have a direction
	bcw = bcw|(onfull&~r);
	bce = bce|(onfull&~l);
	bcs = bcs|(onfull&~u);
	bcn = bcn|(onfull&~d);
	
	%corners--is in two of the previous or is surrounded
	bcc = (bcw&bce)|(bcw&bcs)|(bcw&bcn)|(bce&bcs)|(bce&bcn)|(bcs&bcn);
	
	%inner corner boundary condition
	bcci = onfull&(r&l&u&d);
	
	%bcw = bcw|bcci;
	%bce = bce|bcci;
	%bcs = bcs|bcci;
	%bcn = bcn|bcci;
	bcc = bcc|bcci;
	
	dbcfull.w = bcw;
	dbcfull.e = bce;
	dbcfull.s = bcs;
	dbcfull.n = bcn;
	dbcfull.c = bcc;
	
	%wipe out invalid indices
	if(strcmp(side,'inner'))
		dbc.w = bcw(valindinner);
		dbc.e = bce(valindinner);
		dbc.s = bcs(valindinner);
		dbc.n = bcn(valindinner);
		dbc.c = bcc(valindinner);
	else
		dbc.w = bcw(valindouter);
		dbc.e = bce(valindouter);
		dbc.s = bcs(valindouter);
		dbc.n = bcn(valindouter);
		dbc.c = bcc(valindouter);
	end
	
	
end

