function par = ParametersDebug
	% Parameters
	%
	% Definite parameters go here
	
	addpath('nssolvers');
	addpath('solvers');
	addpath('bcfunctions');
	addpath('rhfunctions');
	addpath('maps');
	addpath('utils');
	
	par.maptype = 'g';
	par.mapfile = 'box.txt';
	par.h = 0.01;
	par.ghostpoints = false;
	par.streamfunction = true;
	par.order = 2;
	par.dt = 1e-2;
	par.tf = 4e-0;
	par.timesteps = round(par.tf/par.dt);
	par.usestagger = true;
	
	
	%flow parameters
	par.inflowAmp = 1;
	par.nu = 1;							%kinematic viscosity
	par.Re = 1e2;							%default value
	
	%plotting parameters
	par.toPlot = 1;						%1==normal, 2==debug
	par.filter = false;
	par.numfilter = 1;
	par.conlines = 30;
	par.zeroout = false;
	par.plot = true;
	
	%domain decomposition parameters
	par.ddrun = false;
	par.ddbounds = {{[0.0,-0.5],[1.5,0.5]},{[1.0,-1.5],[3.5,1.5]},{[3.0,-1.5],[5.0,1.5]}};
	par.ddoverlap = 0.5;
	par.ddmidratio = 0.6;
	par.dditer = 10;
	par.topause = 0;
	
	par.rhfunc = @RHZero;
	par.bcfunc = @BCDrivCav;
	par.solver = @SOBih;
	par.ddsolver = @DDMSch;
	par.nssolver = @NSPrim;
	par.gridmaker = @MakeStaggeredGridsBox;
	
end

