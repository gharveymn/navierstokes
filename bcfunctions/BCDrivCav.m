function [rhs,bcio] = BCDrivCav(grids,filtering,rhs,par,grid)
	
	switch grid
		case 'u'
			[rhs,bcio] = bcu(grids,filtering,rhs,par);
		case 'v'
			[rhs,bcio] = bcv(grids,filtering,rhs,par);
		case 'p'
			[rhs,bcio] = bcp(grids,filtering,rhs,par);
		case 'q'
			[rhs,bcio] = bcq(grids,filtering,rhs,par);
		otherwise
			ME = MException('bcsymchns:invalidParameterException','Invalid value for grid');
			throw(ME)
	end
		
end

function [rhs,bcio] = bcu(grids,filtering,rhs,par)
	
	flowrate = 1;
	
	on = filtering.on;
	ymax = max(grids.ymesh(on));
	
	% driven flow
	drivy = ymax*ones(numel(on),1);
	in = ones(numel(on),flowrate);
	in(~(grids.ymesh==drivy)) = 0;
	
	rhs = rhs + in;
	
	bcio = on&~on;
	
end

function [rhs,bcio] = bcv(grids,filtering,rhs,par)
	on = filtering.on;
	bcio = on&~on;
end

function [rhs,bcio] = bcp(grids,filtering,rhs,par)
	on = filtering.on;
	bcio = on&~on;
end

function [rhs,bcio] = bcq(grids,filtering,rhs,par)
	on = filtering.on;
	bcio = on&~on;
end

