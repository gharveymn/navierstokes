function [clomeshfull,gridsnew,filteringnew,ret] = closure(grids,filtering,side,gp,vect)
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
	% none of this is true, this is a HACK
	% THIS SHOULDNT EVEN WORKKKK
	
	%TODO: specify the boundary with gp, cut down required arguments to h and gp
	%		then use a varargin for mesh, meshfull, CAPSmesh etc.
	%		That should generalize the program for other use
	
	xmeshfull = grids{7};
	ymeshfull = grids{8};
	nx = grids{9};
	ny = grids{10};
	h = grids{11};
	
	valindinner = filtering{2}{1};
	valindouter = filtering{2}{2};
	onfull = filtering{3}{2};
	
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
		
		newmat = zeros(ny+incs+incn,nx+incw+ince);
		newmatl = logical(newmat);
		
		Xmeshnew = reshape(xmeshfull,[nx,ny])';
		Ymeshnew = reshape(ymeshfull,[nx,ny])';
		%Onfull = reshape(onfull,[nx,ny])';
		Valindouter = reshape(valindouter,[nx,ny])';
		Valindinner = reshape(valindinner, [nx,ny])';
		
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
		i2 = i1+ny-1;
		j1 = 1+incw;
		j2 = j1+nx-1;
		
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
		Onfullnew(i1:i2,j1:j2) = reshape(onfull,[nx,ny])';
		Bcw(i1:i2,j1:j2) = reshape(bcw,[nx,ny])';
		Bce(i1:i2,j1:j2) = reshape(bce,[nx,ny])';
		Bcs(i1:i2,j1:j2) = reshape(bcs,[nx,ny])';
		Bcn(i1:i2,j1:j2) = reshape(bcn,[nx,ny])';
		Bcc(i1:i2,j1:j2) = reshape(bcc,[nx,ny])';
		
		Clomeshfull = Clomeshfull | circshift(Bcw,-1,2) | circshift(Bce,1,2) | circshift(Bcs,-1) | circshift(Bcn,1);
		
		Clomeshfull = Clomeshfull | circshift(circshift(Bcc,1),1,2)...
			| circshift(circshift(Bcc,1),-1,2)...
			| circshift(circshift(Bcc,-1),1,2)...
			| circshift(circshift(Bcc,-1),-1,2);
		
		%& ~Valindnew wipes out shifted indices which are inside the polgon
		Clomeshfull = Clomeshfull & ~Valindouternew;
		
		Xmeshnew = Xmeshnew./Valindouternew;
		Ymeshnew = Ymeshnew./Valindouternew;
		
		xinitnew = (xminnew:h:xmaxnew)';
		yinitnew = (yminnew:h:ymaxnew)';
		
		nxnew = numel(xinitnew);
		nynew = numel(yinitnew);
		
		xmeshfullnew = kron(ones(nynew,1),xinitnew);
		ymeshfullnew = kron(yinitnew,ones(nxnew,1));
		
		%clomeshfull is the closure to the polygon
		clomeshfull = reshape(Clomeshfull',[nxnew*nynew,1]);
		
		%onfullnew and its derivatives stay in the same relative position as the grid expands
		onfullnew = reshape(Onfullnew',[nxnew*nynew,1]);
		
		Valindouternew = Valindouternew | Clomeshfull;
		
		%valindnew includes clomeshfull
		valindouternew = reshape(Valindouternew',[nxnew*nynew,1]);
		valindinnernew = reshape(Valindinnernew',[nxnew*nynew,1]);
		
		filterMatnew = spdiag(valindouternew);
		filterMatnew = filterMatnew(valindouternew,:);
		
		xmeshnew = filterMatnew*xmeshfullnew;
		ymeshnew = filterMatnew*ymeshfullnew;
		
		onnew = onfullnew(valindouternew);
		
		if(exist('vect','var'))
			Ret = newmat;
			Ret(i1:i2,j1:j2) = reshape(vect,[nx,ny])';
			ret = reshape(Ret',[nxnew*nynew,1]);
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
		
		newmat = zeros(ny+incs+incn,nx+incw+ince);
		newmatl = logical(newmat);
		
		%Xmeshnew = reshape(xmeshfull,[nx,ny])';
		%Ymeshnew = reshape(ymeshfull,[nx,ny])';
		%Onfull = reshape(onfull,[nx,ny])';
		Valindinner = reshape(valindinner,[nx,ny])';
		Valindouter = reshape(valindouter,[nx,ny])';
		
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
		i2 = i1+ny-1;
		j1 = 1+incw;
		j2 = j1+nx-1;
		
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
		
		Xmeshnew(i1:i2,j1:j2) = reshape(xmeshfull,[nx,ny])';
		Ymeshnew(i1:i2,j1:j2) = reshape(ymeshfull,[nx,ny])';
		Valindinnernew(i1:i2,j1:j2) = Valindinner;
		Valindouternew(i1:i2,j1:j2) = Valindouter;
		Gpmatinds(i1:i2,j1:j2) = reshape(gp,[nx,ny])';
		Onfullnew(i1:i2,j1:j2) = reshape(onfull,[nx,ny])';
		Bcw(i1:i2,j1:j2) = reshape(bcw,[nx,ny])';
		Bce(i1:i2,j1:j2) = reshape(bce,[nx,ny])';
		Bcs(i1:i2,j1:j2) = reshape(bcs,[nx,ny])';
		Bcn(i1:i2,j1:j2) = reshape(bcn,[nx,ny])';
		Bcc(i1:i2,j1:j2) = reshape(bcc,[nx,ny])';
		Origmatinds(i1:i2,j1:j2) = ones(ny,nx);
		
		nxnew = nx - incw - ince;
		nynew = ny - incn - incs;
		
		i1p = i1+incs;
		i2p = i1p+nynew-1;
		j1p = j1+incw;
		j2p = j1p+nxnew-1;
		
		Innermatinds(i1p:i2p,j1p:j2p) = ones(nynew,nxnew);
		
		Clomeshfull = Clomeshfull | circshift(Bcw,1,2) | circshift(Bce,-1,2) | circshift(Bcs,1) | circshift(Bcn,-1);
		
		Clomeshfull = Clomeshfull | circshift(circshift(Bcc,1),1,2)...
			| circshift(circshift(Bcc,1),-1,2)...
			| circshift(circshift(Bcc,-1),1,2)...
			| circshift(circshift(Bcc,-1),-1,2);
		
		Valindinnernew = Valindinnernew & ~Gpmatinds;
		Valindouternew = Valindouternew & ~Gpmatinds;
		
		%& Valindnew wipes out shifted indices which are outside the polgon
		Clomeshfull = Clomeshfull & Valindinnernew;
		
		Xmeshnew = Xmeshnew./Valindinnernew;
		Ymeshnew = Ymeshnew./Valindinnernew;
		
		Clomeshfull = reshape(Clomeshfull(Origmatinds),[ny,nx]);
		Valindinnernew = reshape(Valindinnernew(Innermatinds),[nynew,nxnew]);
		Valindouternew = reshape(Valindouternew(Innermatinds),[nynew,nxnew]);
		Onmeshfull = reshape(Onfullnew(Innermatinds),[nynew,nxnew]);
		Xmeshnew = reshape(Xmeshnew(Innermatinds),[nynew,nxnew]);
		Ymeshnew = reshape(Ymeshnew(Innermatinds),[nynew,nxnew]);
		
		xinitnew = (xminnew-eps:h:xmaxnew+eps)';
		yinitnew = (yminnew-eps:h:ymaxnew+eps)';
		
		xmeshfullnew = kron(ones(nynew,1),xinitnew);
		ymeshfullnew = kron(yinitnew,ones(nxnew,1));
		
		%clomeshfull remains the same size as the original grids --- this is for using the function in selection mode
		clomeshfull = reshape(Clomeshfull',[nx*ny,1]);
		
		%onfullnew is the border on the new smaller domain
		onfullnew = reshape(Onmeshfull',[nxnew*nynew,1]);
		
		Valindinnernew = Valindinnernew | Onmeshfull;
		valindinnernew = reshape(Valindinnernew',[nxnew*nynew,1]);
		valindouternew = reshape(Valindouternew',[nxnew*nynew,1]);
		
		filterMatnew = spdiag(valindinnernew);
		filterMatnew = filterMatnew(valindinnernew,:);
		
		xmeshnew = filterMatnew*xmeshfullnew;
		ymeshnew = filterMatnew*ymeshfullnew;
		
		onnew = onfullnew(valindouternew);
		
		
		if(exist('vect','var'))
			Ret = newmat;
			Ret(i1:i2,j1:j2) = reshape(vect,[nx,ny])';
			Ret = reshape(Ret(Innermatinds),[nynew,nxnew]);
			ret = reshape(Ret',[nxnew*nynew,1]);
		else
			ret = [];
		end
		
	else
		ME = MException('closure:invalidParameterException','Invalid value for side');
		throw(ME)
	end
	
	gridsnew = {xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshfullnew,ymeshfullnew,nxnew,nynew,h};		
	filteringnew = {filterMatnew,{valindinnernew,valindouternew},{onnew,onfullnew},{dbc,dbcfull}};
	[dbc,dbcfull] = boundarysides(gridsnew,filteringnew);
	filteringnew{4} = {dbc,dbcfull};
	
	
	
	
	
end

