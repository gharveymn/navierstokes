function par = Parameters
	% Parameters
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
	par.dt = 0.1;
	par.tf = 200;
	
	%flow parameters
	par.inflowAmp = 1;
	par.nu = 1;							%kinematic viscosity
	par.Re = 1e2;							%default value
	par.omega = 0.5;
	
	%plotting parameters
	par.toPlot = 2;						%1==normal, 2==debug, 3==special
	par.conlines = 40;
	par.quivVectSca = .1*(par.h/0.05);
	par.plotoniter = 100;
	par.noplot = false;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymChNS;
	par.nssolver = @NSPrim;
	par.gridmaker = @MakeStaggeredGrids;
	par.model = @symchlong;

	%applies a specific model
	par = par.model(par);
	
	par.timesteps = round(par.tf/par.dt);
	
	
end

