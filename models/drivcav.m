function par = drivcav(par)
	
	par.mapfile = 'box.txt';
	par.bcfunc = @BCDrivCav;
	par.h = 0.01;
	
end

