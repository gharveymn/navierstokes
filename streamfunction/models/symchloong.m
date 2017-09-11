function par = symchloong(par)
	
	par.mapfile = 'symchloong.txt';
	par.bcfunc = @BCConIO;
	par.hx = 0.25;
	par.hy = 0.1;
	par.Re = 45;
	par.plotoniter = 100;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.1;
	par.tf = 1000;
	
end