function [rhs,bc] = BCZero(grids,filtering,rhs,par)
	bc = {{filtering{3}{1},filtering{3}{1}},{filtering{3}{2},filtering{3}{2}}};
	rhs(bc{1}{1}|bc{1}{2}) = 0;
end