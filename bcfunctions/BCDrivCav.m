function [rhs,bcio] = BCSymChNS(grids,filtering,rhs,par)
	
	flowrate = 1;
	
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	nxp1 = par.nxp1;
	valind = filtering.valind;
	
	on = filtering.on;
	onfull = filtering.onfull;
	dbcfull = filtering.dbcfull;
	gpca = filtering.gp;
	
	h = par.h;
	
	xmax = max(xmesh(on));
	xmin = min(xmesh(on));
	ymax = max(ymesh(on));
	ymin = min(ymesh(on));
	
	% driven flow
	drivy = ymax*ones(numel(on),1);
	in = ones(numel(on),flowrate);
	in(~(ymesh==drivy)) = 0;
	
	rhs = rhs + in;
	
	bcio = on&~on;
end

