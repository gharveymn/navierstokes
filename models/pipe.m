function par = pipe(par)
	par.mapfile = 'rect.txt';
	par.bcfunc = @BCConIO;
	par.hx = 0.05;
	par.hy = 0.01;
	par.Re = 1000;
	par.plotoniter = 10;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.01;
	par.tf = 1000;
	
end
