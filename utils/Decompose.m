function [gridsnew,filteringnew,qmeshnew,bcnew,rhsnew,bcb] = Decompose(grids,filtering,qmesh,rhs,par)
	%DECOMPOSE for use with symch
	
	ddbounds = par.ddbounds;
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	del = (par.order-1)*par.h+par.h/2; %the half h is for for floating point arithmetic errors
	
	[g1,f1,p1] = getGrid(grids,filtering,qmesh,ddbounds{1},par,del);
	[g2,f2,p2] = getGrid(grids,filtering,qmesh,ddbounds{2},par,del);
	[g3,f3,p3] = getGrid(grids,filtering,qmesh,ddbounds{3},par,del);
	
	gridsnew = {g1,g2,g3};
	filteringnew = {f1,f2,f3};
	qmeshnew = {p1,p2,p3};
	
	bc1 = f1{3}{1};
	bc2 = f2{3}{1};
	bc3 = f3{3}{1};
	
	for i=1:par.order-1
		bc1 = bc1 | f1{5}{i};
		bc2 = bc2 | f2{5}{i};
		bc3 = bc3 | f3{5}{i};
	end
	
	
	%bounding lines
	blx11 = ddbounds{1}{1}(1);		%inlet
	blx12 = ddbounds{1}{2}(1);		%b1
	bly11 = ddbounds{1}{1}(2);		%lower y bound
	bly12 = ddbounds{1}{2}(2);		%upper y bound
	
	blx21 = ddbounds{2}{1}(1);		%b0
	blx22 = ddbounds{2}{2}(1);		%b3
	bly21 = ddbounds{2}{1}(2);		%lower y bound
	bly22 = ddbounds{2}{2}(2);		%upper y bound
	
	blx31 = ddbounds{3}{1}(1);		%b2
	blx32 = ddbounds{3}{2}(1);		%outlet
	bly31 = ddbounds{3}{1}(2);		%lower y bound
	bly32 = ddbounds{3}{2}(2);		%upper y bound
	
	%overlaps are b0 <= x <= b1 and b2 <= x < b3
	
	xmesh1 = g1{3};
	ymesh1 = g1{4};
	
	xmesh2 = g2{3};
	ymesh2 = g2{4};
	
	xmesh3 = g3{3};
	ymesh3 = g3{4};
	
	%notation: bcx{y} means boundary x located belonging to the region y
	%this is because we need indices for both in each case
	
	
	%bc0 is all with x <= b1 inside R1 on boundary of R2
	bcb0{1} = (xmesh1 >= blx21-del) & (xmesh1 <= blx21+eps) & (ymesh1 >= bly21-del) & (ymesh1 <= bly22+del) ...
											  & (ymesh1 >= bly11-del) & (ymesh1 <= bly12+del);
	bcb0{2} = bc2 & (xmesh2 <= blx21+eps) & (ymesh2 >= bly11-del) & (ymesh2 <= bly12+del);
	
	%belong to R2
	%bc1 is all with x >= b0 inside R2 on boundary of R1
	bcb1{2} = ((xmesh2 <= blx12+del) & (xmesh2 >= blx12-eps) & (ymesh2 >= bly11-del) & (ymesh2 <= bly12+del)) ...
						| ((xmesh2 <= blx12+del) & (((ymesh2 >= bly11-del) & (ymesh2 <= bly11+eps)) ...
											 | ((ymesh2 >= bly12-eps) & (ymesh2 <= bly12+del))));
	bcb1{1} = bc1 & (xmesh1 >= blx21-del) & (ymesh1 >= bly21-del) & (ymesh1 <= bly22+del);
	
	%bc2 is all with x <= b3 inside R2 on boundary of R3
	bcb2{2} = (xmesh2 >= blx31-del) & (xmesh2 <= blx31+eps) & (ymesh2 >= bly31-del) & (ymesh2 <= bly32+del);
	bcb2{3} = bc3 & (xmesh3 <= blx31+eps) & (ymesh3 >= bly21-del) & (ymesh3 <= bly22+del);
	
	%bc3 is all with x >= b2 inside R3 on boundary of R2
	bcb3{3} = (xmesh3 >= blx22-del) & (xmesh3 <= blx22+eps) & (ymesh3 >= bly21-del) & (ymesh3 <= bly22+del);
	bcb3{2} = bc2 & (xmesh2 >= blx22-eps) & (ymesh2 >= bly21-del) & (ymesh2 <= bly22+del);
	
	bcsur = filtering.on;
	for i=1:par.order-1
		bcsur = bcsur | filtering.gp{i}(filtering.valindouter);
	end
	
	onouter1 = bcsur((xmesh >= blx11-del) & (xmesh <= blx12+del) & (ymesh >= bly11-del) & (ymesh <= bly12+del));
	onouter2 = bcsur((xmesh >= blx21-del) & (xmesh <= blx22+del) & (ymesh >= bly21-del) & (ymesh <= bly22+del));
	onouter3 = bcsur((xmesh >= blx31-del) & (xmesh <= blx32+del) & (ymesh >= bly31-del) & (ymesh <= bly32+del));

	rhs1 = rhs((xmesh >= blx11-del) & (xmesh <= blx12+del) & (ymesh >= bly11-del) & (ymesh <= bly12+del));
	rhs2 = rhs((xmesh >= blx21-del) & (xmesh <= blx22+del) & (ymesh >= bly21-del) & (ymesh <= bly22+del));
	rhs3 = rhs((xmesh >= blx31-del) & (xmesh <= blx32+del) & (ymesh >= bly31-del) & (ymesh <= bly32+del));		
	

	bc1 = bcb1{1} | onouter1;
	bc2 = bcb0{2} | bcb3{2} | onouter2;
	bc3 = bcb2{3} | onouter3;

	%for now
	bc1x = bc1;
	bc1y = bc1;
	bc2x = bc2;
	bc2y = bc2;
	bc3x = bc3;
	bc3y = bc3;
	
	bcnew = {{{bc1x,bc1y},{bc1x,bc1y}},{{bc2x,bc2y},{bc2x,bc2y}},{{bc3x,bc3y},{bc3x,bc3y}}};
	rhsnew = {rhs1,rhs2,rhs3};
	bcb = {bcb0,bcb1,bcb2,bcb3};

end

function [g,f,qmeshnew] = getGrid(grids,filtering,qmesh,bnds,par,del)
	
	xinit = grids.xinit;
	yinit = grids.yinit;
	xmesh = grids.xmesh;
	ymesh = grids.ymesh;
	
	vinds = (xmesh >= bnds{1}(1)-del) & (xmesh <= bnds{2}(1)+del) ...
				& (ymesh >= bnds{1}(2)-del) & (ymesh <= bnds{2}(2)+del);
	
	%minds = (Xmesh <= bnds{1}(1)) & (Xmesh <= bnds{2}(1)) & (Ymesh <= bnds{1}(2)) & (Ymesh <= bnds{2}(2));
	xinitnew = xinit((xinit >= bnds{1}(1)-del)&(xinit <= bnds{2}(1)+del));
	yinitnew = yinit((yinit >= bnds{1}(2)-del)&(yinit <= bnds{2}(2)+del));
	
	nxnew = numel(xinitnew);
	nynew = numel(yinitnew);
	
	xmeshnew = xmesh(vinds);
	ymeshnew = ymesh(vinds);
	qmeshnew = qmesh(vinds);
	
	Xmeshnew = reshape(xmeshnew,[nxnew,nynew])';
	Ymeshnew = reshape(ymeshnew,[nxnew,nynew])';
	
	g = {xinitnew,yinitnew,xmeshnew,ymeshnew,Xmeshnew,Ymeshnew,xmeshnew,ymeshnew,nxnew,nynew,grids.h};

	%filterMat,{valindinner,valindouter},{on,onfull},{dbc,dbcfull},{gp1,gp2}
	
	%these are just rectangles, so we can do this
	filterMatnew = speye(nxnew*nynew);
	valindinnernew = ones(size(xmeshnew));
	valindouternew = valindinnernew;
	onnew = kron([1;zeros(nynew-2,1);1],ones(nxnew,1))|kron(ones(nynew,1),[1;zeros(nxnew-2,1);1]);
	
	f = {filterMatnew,{valindinnernew,valindouternew},{onnew,onnew},{onnew,onnew},{[],[]}};
	%[dbc,dbcfull] = boundarysides(grids,filtering,gp,'outer');

	if(par.ghostpoints)
		
		switch par.order
			case 1
				%do nothing
			case 2
				gp1 = f{3}{1};
				on = closure(g,f,'inner');
				f{5} = {gp1};
				f{3} = {on,on};
			case 3
				[gp1,gtmp,ftmp] = closure(g,f,'inner');
				[gp2,gtmp,ftmp] = closure(gtmp,ftmp,'inner');
				[~,~,~,gp2] = closure(gtmp,ftmp,[],[],gp2,gp2);
				f = [f,{gp2}];
				f{5} = {gp1,gp2};
			otherwise
				ME = MException('closure:invalidParameterException','Invalid value for par.order');
				throw(ME)
		end
		
	end
	
	
	
end

