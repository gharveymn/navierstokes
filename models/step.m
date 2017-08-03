function par = step(par)
	
	par.mapfile = 'step.txt';
	par.bcfunc = @BCConIO;
	par.hx = 0.25;
	par.hy = 0.05;
	par.Re = 100;
	par.plotoniter = 1000;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.01;
	par.tf = 1000;
	
end