function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.hx = 0.05;
	par.hy = 0.05;
	par.dt = 0.017;
	par.tf = 20;
	par.plotoniter = 1;
	par.quivVectSca = 1;
	par.steps = 1000;
	
end

