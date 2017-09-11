<<<<<<< HEAD
function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.hx = 0.1;
	par.hy = 0.1;
	par.dt = 0.017;
	par.tf = 20;
	par.plotoniter = 1;
	par.quivVectSca = 1;
	
end

=======
function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.hx = 0.05;
	par.hy = 0.05;
	par.dt = 0.017;
	par.tf = 20;
	par.plotoniter = 1;
	par.quivVectSca = 1;
	
end

>>>>>>> ad23021a4a184400fc53cc4b3919faa5be8ae99f
