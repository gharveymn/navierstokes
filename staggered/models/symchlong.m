<<<<<<< HEAD
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
	
=======
function par = symchlong(par)
	
	par.mapfile = 'symchlong.txt';
	par.bcfunc = @BCConIO;
	par.hx = 0.25;
	par.hy = 0.1;
	par.pnx = 40;
	par.pny = 30;
	par.Re = 30;
	par.toplot = 4;
	par.plotoniter = 1;
	par.conlines = 40;
	par.quivVectSca = .1*(sqrt(par.hx^2 + par.hy^2)/0.05);
	par.dt = 0.1;
	par.tf = 500;
	par.useinterp = true;
	
>>>>>>> ad23021a4a184400fc53cc4b3919faa5be8ae99f
end