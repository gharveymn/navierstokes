function par = symch(par)
	
	par.mapfile = 'symch.txt';
	par.bcfunc = @BCSymChNS;
	par.h = 0.05;
	par.Re = 1e3;
	par.plotoniter = 100;
	par.quivVectSca = .5;
	
end

