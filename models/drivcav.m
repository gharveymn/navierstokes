function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.h = 0.01;
	par.dt = 0.1;
	par.tf = 20;
	par.plotoniter = 1;
	par.quivVectSca = .5;
	
end

