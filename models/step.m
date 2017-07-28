function par = step(par)
	
	par.mapfile = 'step.txt';
	par.bcfunc = @BCSymChNS;
	par.h = 0.05;
	
end