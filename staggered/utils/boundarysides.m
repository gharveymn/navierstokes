function [dbc,dbcfull] = boundarysides(grids,filtering,par,side,nx)
	%BOUNDARYSIDES finds which side each boundary is on
	
	if(~exist('side','var'))
		side = 'inner';
	end
	
	xmeshfull = grids.(side).xmeshfull;
	ymeshfull = grids.(side).ymeshfull;
	onfull = filtering.(side).onfull;
	valind = filtering.(side).valind;
	
	wesnc = [par.wesn {'c','co','ci'}];
	
	xmin = min(xmeshfull);
	xmax = max(xmeshfull);
	ymin = min(ymeshfull);
	ymax = max(ymeshfull);
	
	%get indices on the ends
	xminb = (xmeshfull==xmin);
	xmaxb = (xmeshfull==xmax);
	yminb = (ymeshfull==ymin);
	ymaxb = (ymeshfull==ymax);
	
	%if an index is on the overall boundary then it must be a boundary point
	bcw = onfull&xminb;
	bce = onfull&xmaxb;
	bcs = onfull&yminb;
	bcn = onfull&ymaxb;
	
	%make shifted index matrices, filter out those indices at the max values since circshift loops
	r = circshift(valind&~xmaxb,1);
	l = circshift(valind&~xminb,-1);
	u = circshift(valind&~ymaxb,nx);
	d = circshift(valind&~yminb,-nx);
	
	%if the shift makes it go off the boundary then we have a direction
	bcw = bcw|(onfull&~r);
	bce = bce|(onfull&~l);
	bcs = bcs|(onfull&~u);
	bcn = bcn|(onfull&~d);
	
	%corners--is in two of the previous or is surrounded
	bcco = (bcw&bce)|(bcw&bcs)|(bcw&bcn)|(bce&bcs)|(bce&bcn)|(bcs&bcn);
	
	%inner corner boundary condition
	bcci = onfull&(r&l&u&d);
	
	%bcw = bcw|(bcci&(circshift(bcw,-nx)|circshift(bcw,nx)));
	%bce = bce|(bcci&(circshift(bce,-nx)|circshift(bce,nx)));
	%bcs = bcs|(bcci&(circshift(bcs,-1)|circshift(bcs,1)));
	%bcn = bcn|(bcci&(circshift(bcn,-1)|circshift(bcn,1)));
	bcc = bcco|bcci;
	
	
	for i = 1:numel(wesnc)
		eval(['dbcfull.' wesnc{i} '= bc' wesnc{i} ';']);
		
		%wipe out invalid indices
		eval(['dbc.' wesnc{i} '= bc' wesnc{i} '(valind);']);
	end
	
	
end

