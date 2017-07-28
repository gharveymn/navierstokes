function par = symchlong(par)
	
	par.mapfile = 'symchlong.txt';
	par.bcfunc = @BCSymChNS;
	par.hx = 0.25;
	par.hy = 0.25;
	par.Re = 45;
	par.plotoniter = 100;
	par.quivVectSca = .1*(par.h/0.05);
	par.dt = 0.1;
	par.tf = 200;
	
end