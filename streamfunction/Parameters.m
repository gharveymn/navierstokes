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
	
	%set true to switch to ParametersDebug
	par.useGPU = false;
	
	par.varnames = {'u','v','p','q'};
	par.wesn = {'w','e','s','n'};
	
	par.maptype = 'g';
	par.mapfile = 'symch.txt';
	par.h = 0.05;
	par.hx = 0.05;
	par.hy = 0.05;
	par.nx = 100;
	par.ny = 30;
	par.griduse = 'h';
	par.streamfunction = true;
	par.order = 2;
	par.usestagger = true;
	par.steps = 1000;
	
	%flow parameters
	par.inflowAmp = 1;
	par.nu = 1;							%kinematic viscosity
	par.Re = 1e2;							%default value
	par.omega = 0.1;
	
	%plotting parameters
	par.toPlot = 2;						%1==normal, 2==debug, 3==special
	par.filter = false;
	par.numfilter = 1;
	par.conlines = 40;
	par.zeroout = false;
	par.quivVectSca = .1*(par.h/0.05);
	par.plotoniter = 100;
	par.noplot = false;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymChNS;
	par.nssolver = @NSIter;
	par.gridmaker = @MakeGrids;
	par.model = @drivcav;
	

	par = par.model(par);
	
	
end

