function [clomeshfull,gridsnew,filteringnew,ret] = closure(grids,filtering,side,gp,vect,leaveCornersAlone)
	%CLOSURE surrounds the input grid with a closure
	% expects vector containing x coordinates and vector containing y coordinates of form
	% 11,21,31,41,12,22,32,42,etc.
	% ie. x is the primary iterator
	%
	% note: expects both x and y to iterate increasing
	%
	% args:
	%	grids		the grids cell array
	%	filtering		the filtering cell array
	%	side			specifies whether to surround the polygon or to get one grid in from the boundary
	%	gp			specifies the boundary we want to use (since <code>on</code> might not be the actual boundary)
	%	vect			specifies a vector to transform (must be a meshfull)
	%
	%return:
	%	clomeshfull	the closure (note that this is not added to the grid)
	%	gridsnew		{xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshfullnew,ymeshfullnew}
	%	filteringnew	filtering matrices
	%	ret			vect adjusted to the new domain
	%
	
	%TODO: specify the boundary with gp, cut down required arguments to h and gp
	%		then use a varargin for mesh, meshfull, CAPSmesh etc.
	%		That should generalize the program for other use
	
	xmeshfull = grids.inner.xmeshfull;
	ymeshfull = grids.inner.ymeshfull;
	nxp1 = grids.nxp1;
	nyp1 = grids.nyp1;
	h = grids.h;
	
	valindinner = filtering.valindinner;
	valindouter = filtering.valindouter;
	onfull = filtering.onfull;
	
	if(~exist('side','var') || isempty(side))
		side = 'outer';
	end
	
	if(~exist('gp','var') || isempty(gp))
		gp = onfull;
	end
	
	xmax = max(xmeshfull);
	xmin = min(xmeshfull);
	ymax = max(ymeshfull);
	ymin = min(ymeshfull);
	
	if(strcmp(side,'outer'))
		
		[dbc,dbcfull] = boundarysides(grids,filtering,gp,'outer');
		bcw = dbcfull{1};
		bce = dbcfull{2};
		bcs = dbcfull{3};
		bcn = dbcfull{4};
		bcc = dbcfull{5};
		
		% how much to increase the grid size in NSEW directions
		incw = ~isempty(xmeshfull(xmeshfull==xmin & valindouter));
		ince = ~isempty(xmeshfull(xmeshfull==xmax & valindouter));
		incs = ~isempty(ymeshfull(ymeshfull==ymin & valindouter));
		incn = ~isempty(ymeshfull(ymeshfull==ymax & valindouter));
		
		newmat = zeros(nyp1+incs+incn,nxp1+incw+ince);
		newmatl = logical(newmat);
		
		Xmeshnew = reshape(xmeshfull,[nxp1,nyp1])';
		Ymeshnew = reshape(ymeshfull,[nxp1,nyp1])';
		%Onfull = reshape(onfull,[nx,ny])';
		Valindouter = reshape(valindouter,[nxp1,nyp1])';
		Valindinner = reshape(valindinner, [nxp1,nyp1])';
		
		%NOTE: "SOUTH" INDICES ARE ACTUALLY AT THE TOP OF THE MATRIX
		
		xminnew = xmin;
		xmaxnew = xmax;
		yminnew = ymin;
		ymaxnew = ymax;
		
		if(incw)
			xminnew = xmin - h;
			Xmeshnew = horzcat(xminnew*ones(size(Xmeshnew,1),1),Xmeshnew);
			Ymeshnew = horzcat(Ymeshnew(:,1),Ymeshnew);
		end
		
		if(incn)
			xmaxnew = xmax + h;
			Xmeshnew = horzcat(Xmeshnew,xmaxnew*ones(size(Xmeshnew,1),1));
			Ymeshnew = horzcat(Ymeshnew,Ymeshnew(:,end));
		end
		
		if(incs)
			yminnew = ymin - h;
			Xmeshnew = vertcat(Xmeshnew(1,:),Xmeshnew);
			Ymeshnew = vertcat(yminnew*ones(1,size(Ymeshnew,2)),Ymeshnew);
		end
		
		if(incn)
			ymaxnew = ymax + h;
			Xmeshnew = vertcat(Xmeshnew,Xmeshnew(end,:));
			Ymeshnew = vertcat(Ymeshnew,ymaxnew*ones(1,size(Ymeshnew,2)));
		end
		
		i1 = 1+incs;
		i2 = i1+nyp1-1;
		j1 = 1+incw;
		j2 = j1+nxp1-1;
		
		%Xmeshnew = newmat;
		%Ymeshnew = newmat;
		Clomeshfull = newmatl;
		Onfullnew = newmatl;
		Valindouternew = newmatl;
		Valindinnernew = newmatl;
		Bcw = newmatl;
		Bce = newmatl;
		Bcs = newmatl;
		Bcn = newmatl;
		Bcc = newmatl;
		
		%Xmeshnew(i1:i2,j1:j2) = reshape(xmeshfull,[nx,ny])';
		%Ymeshnew(i1:i2,j1:j2) = reshape(ymeshfull,[nx,ny])';
		Valindouternew(i1:i2,j1:j2) = Valindouter;
		Valindinnernew(i1:i2,j1:j2) = Valindinner;
		Onfullnew(i1:i2,j1:j2) = reshape(onfull,[nxp1,nyp1])';
		Bcw(i1:i2,j1:j2) = reshape(bcw,[nxp1,nyp1])';
		Bce(i1:i2,j1:j2) = reshape(bce,[nxp1,nyp1])';
		Bcs(i1:i2,j1:j2) = reshape(bcs,[nxp1,nyp1])';
		Bcn(i1:i2,j1:j2) = reshape(bcn,[nxp1,nyp1])';
		Bcc(i1:i2,j1:j2) = reshape(bcc,[nxp1,nyp1])';
		
		Clomeshfull = Clomeshfull | circshift(Bcw,-1,2) | circshift(Bce,1,2) | circshift(Bcs,-1) | circshift(Bcn,1);
		
		if(~exist('leaveCornersAlone','var') || ~leaveCornersAlone)
			Clomeshfull = Clomeshfull | circshift(circshift(Bcc,1),1,2)...
				| circshift(circshift(Bcc,1),-1,2)...
				| circshift(circshift(Bcc,-1),1,2)...
				| circshift(circshift(Bcc,-1),-1,2);
		end
		
		%& ~Valindnew wipes out shifted indices which are inside the polgon
		Clomeshfull = Clomeshfull & ~Valindouternew;
		
		Xmeshnew = Xmeshnew./Valindouternew;
		Ymeshnew = Ymeshnew./Valindouternew;
		
		xinitnew = (xminnew:h:xmaxnew)';
		yinitnew = (yminnew:h:ymaxnew)';
		
		nxp1new = numel(xinitnew);
		nyp1new = numel(yinitnew);
		
		xmeshfullnew = kron(ones(nyp1new,1),xinitnew);
		ymeshfullnew = kron(yinitnew,ones(nxp1new,1));
		
		%clomeshfull is the closure to the polygon
		clomeshfull = reshape(Clomeshfull',[nxp1new*nyp1new,1]);
		
		%onfullnew and its derivatives stay in the same relative position as the grid expands
		onfullnew = reshape(Onfullnew',[nxp1new*nyp1new,1]);
		
		Valindouternew = Valindouternew | Clomeshfull;
		
		%valindnew includes clomeshfull
		valindouternew = reshape(Valindouternew',[nxp1new*nyp1new,1]);
		valindinnernew = reshape(Valindinnernew',[nxp1new*nyp1new,1]);
		
		filterMatnew = spdiag(valindouternew);
		filterMatnew = filterMatnew(valindouternew,:);
		
		xmeshnew = filterMatnew*xmeshfullnew;
		ymeshnew = filterMatnew*ymeshfullnew;
		
		onnew = onfullnew(valindouternew);
		
		if(exist('vect','var'))
			Ret = newmat;
			Ret(i1:i2,j1:j2) = reshape(vect,[nxp1,nyp1])';
			ret = reshape(Ret',[nxp1new*nyp1new,1]);
		else
			ret = [];
		end
		
	elseif(strcmp(side,'inner'))
		%returns smaller size
		
		%clomeshfull should be wrt the original mesh
		%everything else should be converted to the smaller size
		
		[dbc,dbcfull] = boundarysides(grids,filtering,gp,'outer');
		bcw = dbcfull{1};
		bce = dbcfull{2};
		bcs = dbcfull{3};
		bcn = dbcfull{4};
		bcc = dbcfull{5};
		
		% how much to increase the grid size in NSEW directions
		incw = ~isempty(xmeshfull(xmeshfull==xmin & valindouter));
		ince = ~isempty(xmeshfull(xmeshfull==xmax & valindouter));
		incs = ~isempty(ymeshfull(ymeshfull==ymin & valindouter));
		incn = ~isempty(ymeshfull(ymeshfull==ymax & valindouter));
		
		newmat = zeros(nyp1+incs+incn,nxp1+incw+ince);
		newmatl = logical(newmat);
		
		%Xmeshnew = reshape(xmeshfull,[nx,ny])';
		%Ymeshnew = reshape(ymeshfull,[nx,ny])';
		%Onfull = reshape(onfull,[nx,ny])';
		Valindinner = reshape(valindinner,[nxp1,nyp1])';
		Valindouter = reshape(valindouter,[nxp1,nyp1])';
		
		%NOTE: "SOUTH" INDICES ARE ACTUALLY AT THE TOP OF THE MATRIX
		
		xminnew = xmin;
		xmaxnew = xmax;
		yminnew = ymin;
		ymaxnew = ymax;
		
		if(incw)
			xminnew = xmin + h;
			%Xmeshnew = horzcat(xminnew*ones(ny,1),Xmeshnew);
			%Ymeshnew = horzcat(Ymeshnew(:,1),Ymeshnew);
		end
		
		if(incn)
			xmaxnew = xmax - h;
			%Xmeshnew = horzcat(Xmeshnew,xmaxnew*ones(ny,1));
			%Ymeshnew = horzcat(Ymeshnew,Ymeshnew(:,end));
		end
		
		if(incs)
			yminnew = ymin + h;
			%Xmeshnew = vertcat(Xmeshnew(1,:),Xmeshnew);
			%Ymeshnew = vertcat(yminnew*ones(1,nx),Ymeshnew);
		end
		
		if(incn)
			ymaxnew = ymax - h;
			%Xmeshnew = vertcat(Xmeshnew,Xmeshnew(end,:));
			%Ymeshnew = vertcat(Ymeshnew,ymaxnew*ones(1,nx));
		end
		
		i1 = 1+incs;
		i2 = i1+nyp1-1;
		j1 = 1+incw;
		j2 = j1+nxp1-1;
		
		Xmeshnew = newmat;
		Ymeshnew = newmat;
		Clomeshfull = newmatl;
		Valindinnernew = newmatl;
		Valindouternew = newmatl;
		Onfullnew = newmatl;
		Bcw = newmatl;
		Bce = newmatl;
		Bcs = newmatl;
		Bcn = newmatl;
		Bcc = newmatl;
		Origmatinds = newmatl;
		Innermatinds = newmatl;
		Gpmatinds = newmatl;
		
		Xmeshnew(i1:i2,j1:j2) = reshape(xmeshfull,[nxp1,nyp1])';
		Ymeshnew(i1:i2,j1:j2) = reshape(ymeshfull,[nxp1,nyp1])';
		Valindinnernew(i1:i2,j1:j2) = Valindinner;
		Valindouternew(i1:i2,j1:j2) = Valindouter;
		Gpmatinds(i1:i2,j1:j2) = reshape(gp,[nxp1,nyp1])';
		Onfullnew(i1:i2,j1:j2) = reshape(onfull,[nxp1,nyp1])';
		Bcw(i1:i2,j1:j2) = reshape(bcw,[nxp1,nyp1])';
		Bce(i1:i2,j1:j2) = reshape(bce,[nxp1,nyp1])';
		Bcs(i1:i2,j1:j2) = reshape(bcs,[nxp1,nyp1])';
		Bcn(i1:i2,j1:j2) = reshape(bcn,[nxp1,nyp1])';
		Bcc(i1:i2,j1:j2) = reshape(bcc,[nxp1,nyp1])';
		Origmatinds(i1:i2,j1:j2) = ones(nyp1,nxp1);
		
		nxp1new = nxp1 - incw - ince;
		nyp1new = nyp1 - incn - incs;
		
		i1p = i1+incs;
		i2p = i1p+nyp1new-1;
		j1p = j1+incw;
		j2p = j1p+nxp1new-1;
		
		Innermatinds(i1p:i2p,j1p:j2p) = ones(nyp1new,nxp1new);
		
		Clomeshfull = Clomeshfull | circshift(Bcw,1,2) | circshift(Bce,-1,2) | circshift(Bcs,1) | circshift(Bcn,-1);
		
		if(~exist('leaveCornersAlone','var') || ~leaveCornersAlone)
			Clomeshfull = Clomeshfull | circshift(circshift(Bcc,1),1,2)...
				| circshift(circshift(Bcc,1),-1,2)...
				| circshift(circshift(Bcc,-1),1,2)...
				| circshift(circshift(Bcc,-1),-1,2);
		end
		
		Valindinnernew = Valindinnernew & ~Gpmatinds;
		Valindouternew = Valindouternew & ~Gpmatinds;
		
		%& Valindnew wipes out shifted indices which are outside the polgon
		Clomeshfull = Clomeshfull & Valindinnernew;
		
		Xmeshnew = Xmeshnew./Valindinnernew;
		Ymeshnew = Ymeshnew./Valindinnernew;
		
		Clomeshfull = reshape(Clomeshfull(Origmatinds),[nyp1,nxp1]);
		Valindinnernew = reshape(Valindinnernew(Innermatinds),[nyp1new,nxp1new]);
		Valindouternew = reshape(Valindouternew(Innermatinds),[nyp1new,nxp1new]);
		Onmeshfull = reshape(Onfullnew(Innermatinds),[nyp1new,nxp1new]);
		Xmeshnew = reshape(Xmeshnew(Innermatinds),[nyp1new,nxp1new]);
		Ymeshnew = reshape(Ymeshnew(Innermatinds),[nyp1new,nxp1new]);
		
		xinitnew = (xminnew-eps:h:xmaxnew+eps)';
		yinitnew = (yminnew-eps:h:ymaxnew+eps)';
		
		xmeshfullnew = kron(ones(nyp1new,1),xinitnew);
		ymeshfullnew = kron(yinitnew,ones(nxp1new,1));
		
		%clomeshfull remains the same size as the original grids --- this is for using the function in selection mode
		clomeshfull = reshape(Clomeshfull',[nxp1*nyp1,1]);
		
		%onfullnew is the border on the new smaller domain
		onfullnew = reshape(Onmeshfull',[nxp1new*nyp1new,1]);
		
		Valindinnernew = Valindinnernew | Onmeshfull;
		valindinnernew = reshape(Valindinnernew',[nxp1new*nyp1new,1]);
		valindouternew = reshape(Valindouternew',[nxp1new*nyp1new,1]);
		
		filterMatnew = spdiag(valindinnernew);
		filterMatnew = filterMatnew(valindinnernew,:);
		
		xmeshnew = filterMatnew*xmeshfullnew;
		ymeshnew = filterMatnew*ymeshfullnew;
		
		onnew = onfullnew(valindouternew);
		
		
		if(exist('vect','var'))
			Ret = newmat;
			Ret(i1:i2,j1:j2) = reshape(vect,[nxp1,nyp1])';
			Ret = reshape(Ret(Innermatinds),[nyp1new,nxp1new]);
			ret = reshape(Ret',[nxp1new*nyp1new,1]);
		else
			ret = [];
		end
		
	else
		ME = MException('closure:invalidParameterException','Invalid value for side');
		throw(ME)
	end
	
	gridsnew.xinit = xinitnew;
	gridsnew.yinit = yinitnew;
	gridsnew.xmesh = xmeshnew;
	gridsnew.ymesh = ymeshnew;
	gridsnew.Xmesh = Xmeshnew;
	gridsnew.Ymesh = Ymeshnew;
	gridsnew.xmeshfull = xmeshfullnew;
	gridsnew.ymeshfull = ymeshfullnew;
	gridsnew.nxp1 = nxp1new;
	gridsnew.nyp1 = nyp1new;
	gridsnew.h = h;
	
	filteringnew.filterMat = filterMatnew;
	filteringnew.valindinner = valindinnernew;
	filteringnew.valindouter = valindouternew;
	filteringnew.on = onnew;
	filteringnew.onfull = onfullnew;
	filteringnew.dbc = dbc;
	filteringnew.dbcfull = dbcfull;
	
	[dbc,dbcfull] = boundarysides(gridsnew,filteringnew);
	filteringnew.dbc = dbc;
	filteringnew.dbcfull = dbcfull;
	
	
	
	
	
end

