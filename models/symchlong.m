function par = symchlong(par)
	
	par.mapfile = 'symchlong.txt';
	par.bcfunc = @BCSymChNS;
	par.hx = 0.25;
	par.hy = 0.1;
	par.Re = 100;
	par.plotoniter = 10;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.1;
	par.tf = 1000;
	
end