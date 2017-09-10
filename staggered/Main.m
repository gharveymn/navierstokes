if(~exist('testrun','var') || ~testrun)
	par = Parameters;
end

[grids,filtering,par] = MakeStaggeredGrids(par);
[rhs] = par.rhfunc(grids,par);
[rhs,filtering] = par.bcfunc(grids,filtering,rhs,par);

[grids,filtering,res,par] = par.nssolver(par,grids,filtering,rhs);



