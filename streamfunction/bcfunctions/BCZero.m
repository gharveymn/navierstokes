function [rhs,bc] = BCZero(grids,filtering,rhs,par)
	bc = {{filtering.on,filtering.on},{filtering.onfull,filtering.onfull}};
	rhs(bc{1}{1}|bc{1}{2}) = 0;
end