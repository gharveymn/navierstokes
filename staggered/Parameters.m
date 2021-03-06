<<<<<<< HEAD
function par = Parameters
	% PARAMETERS
	%
	% Definite parameters go here
	
	addpath('nssolvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	addpath('utils');
	addpath('models');
	
	%default parameters
	par.useGPU = false;
	
	par.varnames = {'u','v','p','q'};
	par.wesn = {'w','e','s','n'};
	
	par.mapfile = 'symch.txt';
	par.h = 0.05;
	par.hx = 0.05;
	par.hy = 0.05;
	par.nx = 100;
	par.ny = 30;
	par.griduse = 'h';
	par.dt = 0.1;
	par.tf = 200;
	
	
	%flow parameters
	par.inflowAmp = 1;
	par.nu = 1;							%kinematic viscosity
	par.Re = 1e2;							%default value
	par.omega = 0.5;
	
	%plotting parameters
	par.toplot = 2;						%1==normal, 2==surf, 3==contour, %4==simple normal, %5==quiver/contour
	par.plotbc = false;
	par.conlines = 40;
	par.quivVectSca = .1*(par.h/0.05);
	par.plotoniter = 100;
	par.noplot = false;
	par.pnx = 60;
	par.pny = 30;
	par.useinterp = true;
	par.numfps = 15;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymChNS;
	par.nssolver = @NSPrim;
	par.model = @symchlong;

	%applies a specific model (to more easily switch between different parameter sets)
	par = par.model(par);
	
	par.timesteps = round(par.tf/par.dt);
	
	
end

=======
function par = Parameters
	% PARAMETERS
	%
	% Definite parameters go here
	
	addpath('nssolvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	addpath('utils');
	addpath('models');
	
	%default parameters
	par.useGPU = false;
	
	par.varnames = {'u','v','p','q'};
	par.wesn = {'w','e','s','n'};
	
	par.mapfile = 'symch.txt';
	par.h = 0.05;
	par.hx = 0.05;
	par.hy = 0.05;
	par.nx = 100;
	par.ny = 30;
	par.griduse = 'h';
	par.dt = 0.1;
	par.tf = 200;
	
	
	%flow parameters
	par.inflowAmp = 1;
	par.nu = 1;							%kinematic viscosity
	par.Re = 1e2;							%default value
	par.omega = 0.5;
	
	%plotting parameters
	par.toplot = 1;						%1==normal, 2==debug, 3==special, %4==simple
	par.conlines = 40;
	par.quivVectSca = .1*(par.h/0.05);
	par.plotoniter = 100;
	par.noplot = false;
	par.pnx = 60;
	par.pny = 30;
	par.useinterp = true;
	par.numfps = 15;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymChNS;
	par.nssolver = @NSPrim;
	par.model = @symchlong;

	%applies a specific model (to more easily switch between different parameter sets)
	par = par.model(par);
	
	par.timesteps = round(par.tf/par.dt);
	
	
end

>>>>>>> ad23021a4a184400fc53cc4b3919faa5be8ae99f
