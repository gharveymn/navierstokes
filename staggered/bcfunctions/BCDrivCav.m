function [rhs,filtering] = BCDrivCav(grids,filtering,rhs,par,spec)
	if(~exist('spec','var'))
		for i=1:numel(par.varnames)
			
			bcf = str2func(['bc' par.varnames{i}]);
			
			[rhs.(par.varnames{i}).inner,filtering] = bcf(grids,filtering,rhs.(par.varnames{i}).inner,par,'inner');
			[rhs.(par.varnames{i}).outer,filtering] = bcf(grids,filtering,rhs.(par.varnames{i}).outer,par,'outer');
			
		end
	else
		bcf = str2func(['bc' spec.varname]);
		[rhs.(spec.varname).(spec.side),filtering] = bcf(grids,filtering,rhs.(spec.varname).(spec.side),par,spec.side);
	end
	
end

function [rhs,filtering] = bcu(grids,filtering,rhs,par,side)
	
	flowrate = 1;
	
	on = filtering.u.(side).on;
	
	ymax = max(grids.u.(side).yinit);
	
	% driven flow
	drivy = ymax*ones(numel(on),1);
	in = ones(numel(on),flowrate);
	in(~(grids.u.(side).ymesh==drivy)) = 0;
	
	rhs = rhs + in;
	
	filtering.u.(side).bcio = on&~on;
	filtering = makebciofull(filtering,'u',side);
	
end

function [rhs,filtering] = bcv(grids,filtering,rhs,par,side)
	on = filtering.v.(side).on;
	filtering.v.(side).bcio = on&~on;
	filtering = makebciofull(filtering,'v',side);
end

function [rhs,filtering] = bcp(grids,filtering,rhs,par,side)
	on = filtering.p.(side).on;
	filtering.p.(side).bcio = on&~on;
	filtering = makebciofull(filtering,'p',side);
end

function [rhs,filtering] = bcq(grids,filtering,rhs,par,side)
	on = filtering.q.(side).on;
	filtering.q.(side).bcio = on&~on;
	filtering = makebciofull(filtering,'q',side);
end

function [filtering] = makebciofull(filtering,varname,side)
	filtering.(varname).(side).bciofull = logical(filtering.(varname).(side).filterMat'*(1*filtering.(varname).(side).bcio));
end

