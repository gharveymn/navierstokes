function par = symchlong(par)
	
	par.mapfile = 'symchlong.txt';
	par.bcfunc = @BCConIO;
	par.hx = 0.25;
	par.hy = 0.1;
	par.pnx = 40;
	par.pny = 30;
	par.Re = 30;
	par.toplot = 4;
	par.plotoniter = 5;
	par.conlines = 40;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.1;
	par.tf = 500;
	par.useinterp = true;
	
end