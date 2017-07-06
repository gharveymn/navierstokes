function [qmesh,mats] = SOBih(grids,filtering,rhs,bc,mats)
	
	if(nargin == 7)
		M = mats{1};
	else
		nx = grids.nxp1;
		ny = grids.nyp1;
		h = grids.h;
		filterMat = filtering.filterMat;
		
		%make derivative matrices
		bih = biharmonic2(nx,ny,h,bc{1}{2}{1},bc{1}{2}{2});
		M = filterMat*bih*filterMat';
		
		mats = {M};
	end
	
	%disp(['lower bound for condition number: ' num2str(condest(bih))])
	
	qmesh = M\rhs;
	
end

