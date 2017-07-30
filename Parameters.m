function par = Parameters
	% Parameters
	%
	% Definite parameters go here
	
	addpath('nssolvers');
	addpath('solvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	addpath('utils');
	addpath('models');
	
	%default parameters
	
	%set true to switch to ParametersDebug
	par.debug = false;
	par.useGPU = false;
	
	par.varnames = {'u','v','p','q'};
	par.wesn = {'w','e','s','n'};
	
	par.maptype = 'g';
	par.mapfile = 'symch.txt';
	par.h = 0.05;
	par.hx = 0.05;
	par.hy = 0.05;
	par.ghostpoints = false;
	par.streamfunction = true;
	par.order = 2;
	par.dt = 0.01;
	par.tf = 200;
	par.usestagger = true;
	
	%flow parameters
	par.inflowAmp = 1;
	par.nu = 1;							%kinematic viscosity
	par.Re = 1e2;							%default value
	par.omega = 0.5;
	
	%plotting parameters
	par.toPlot = 1;						%1==normal, 2==debug
	par.filter = false;
	par.numfilter = 1;
	par.conlines = 30;
	par.zeroout = false;
	par.plot = true;
	par.quivVectSca = .1*(par.h/0.05);
	par.plotoniter = 100;
	par.noplot = false;
	
	
	%domain decomposition parameters
	par.ddrun = false;
	par.ddbounds = {{[0.0,-0.5],[1.5,0.5]},{[1.0,-1.5],[3.5,1.5]},{[3.0,-1.5],[5.0,1.5]}};
	par.ddoverlap = 0.5;
	par.ddmidratio = 0.6;
	par.dditer = 10;
	par.topause = 0;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCSymChNS;
	par.solver = @SOBih;
	par.ddsolver = @DDMSch;
	par.nssolver = @NSPrim;
	par.gridmaker = @MakeStaggeredGrids;
	par.model = @symchlong;
	
	
	if(par.debug)
		par = ParametersDebug;
	end

	par = par.model(par);
	
	par.timesteps = round(par.tf/par.dt);
	
	
end

