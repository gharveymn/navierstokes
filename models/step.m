function par = step(par)
	
	par.mapfile = 'step.txt';
	par.bcfunc = @BCConIO;
	par.hx = 0.25;
	par.hy = 0.05;
	par.Re = 1000;
	par.plotoniter = 10;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.1;
	par.tf = 1000;
	
end