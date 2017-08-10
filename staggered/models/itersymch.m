function par = itersymch(par)
	
	par.mapfile = 'symch.txt';
	par.bcfunc = @BCSymCh;
	par.gridmaker = @MakeGrids;
	par.nssolver = @NSIter;
	par.usestagger = false;
	par.h = 0.05;
	par.omega = 0.9;
	par.toPlot = 2;
	
end

