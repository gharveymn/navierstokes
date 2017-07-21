function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.h = 0.01;
	par.quivVectSca = .1*(par.h/0.05);
	
end

