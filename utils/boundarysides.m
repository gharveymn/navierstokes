function [dbc,dbcfull] = boundarysides(grids,filtering,gp,side)
	%GETWHEREBOUNDARIES I'm somewhat suprised this actually works
	
	if(~exist('side','var'))
		side = 'inner';
	end
	
	xmeshfull = grids.xmeshfull;
	ymeshfull = grids.ymeshfull;
	
	if(strcmp(side,'outer'))
		onfull = gp;
	else
		onfull = filtering.onfull;
	end
	valindinner = filtering.valindinner;
	valindouter = filtering.valindouter;
	nx = grids.nx;
	
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
	
	if(strcmp(side,'outer'))
		r = circshift(valindouter&~xmaxb,1);
		l = circshift(valindouter&~xminb,-1);
		u = circshift(valindouter&~ymaxb,nx);
		d = circshift(valindouter&~yminb,-nx);
	elseif(strcmp(side,'inner'))
		r = circshift(valindinner&~xmaxb,1);
		l = circshift(valindinner&~xminb,-1);
		u = circshift(valindinner&~ymaxb,nx);
		d = circshift(valindinner&~yminb,-nx);
	else
		ME = MException('boundarysides:invalidParameterException','Invalid value for side');
		throw(ME)
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
	
	dbcfull = {bcw,bce,bcs,bcn,bcc};
	
	%wipe out invalid indices
	bcw = bcw(valindouter);
	bce = bce(valindouter);
	bcs = bcs(valindouter);
	bcn = bcn(valindouter);
	bcc = bcc(valindouter);
	
	dbc = {bcw,bce,bcs,bcn,bcc};
end

