function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.hx = 0.05;
	par.hy = 0.05;
	par.omega = 0.01;
	par.Re = 1e1;
	par.dt = 0.017;
	par.tf = 20;
	par.plotoniter = 100;
	par.quivVectSca = 1;
	par.steps = 10000;
	
end

