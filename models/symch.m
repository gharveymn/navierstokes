function par = symch(par)
	
	par.mapfile = 'symch.txt';
	par.bcfunc = @BCSymChNS;
	par.h = 0.05;
	par.quivVectSca = .1*(par.h/0.05);
	
end

