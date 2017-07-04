if(~exist('par','var'))
	par = Parameters;
end

if(~exist('figs','var'))
	[figs,mat,vec,grids,filtering,bc,rhs] = GetStokesSol(par);
else
	[figs,mat,vec,grids,filtering,bc,rhs] = GetStokesSol(par,figs);
end

[par,mat,vec] = par.nssolver(par,figs,mat,vec,grids,filtering,bc,rhs);

X = mat(:,:,1);
Y = mat(:,:,2);
U = mat(:,:,3);
V = mat(:,:,4);
Q = mat(:,:,5);

x = vec(:,1);
y = vec(:,2);
u = vec(:,3);
v = vec(:,4);
q = vec(:,5);

clear mat vec

if(~exist('testrun','var') || ~testrun)
	clear par
end



