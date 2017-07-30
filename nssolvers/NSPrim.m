function [grids,filtering,res,par,figs] = NSPrim(par,grids,filtering,rhs,figs)
	%NSPRIM primitive formulation of NS
	
	if(~exist('par','var'))
		Main();
		return;
	end
	
	nx = grids.p.nx;
	ny = grids.p.ny;
	hx = par.hx;
	hy = par.hy;
	
	U = zeros([ny,nx-1]);
	V = zeros([ny-1,nx]);
	Ue = zeros([ny+2,nx+1]);
	Ve = zeros([ny+1,nx+2]);
	
	if(par.useGPU)
		U = gpuArray(U);
		V = gpuArray(V);
		Ue = gpuArray(Ue);
		Ve = gpuArray(Ve);
	end
	
	%inner
	udbcfull = filtering.u.inner.dbcfull;
	vdbcfull = filtering.v.inner.dbcfull;
	pdbcfull = filtering.p.inner.dbcfull;
	qdbcfull = filtering.q.inner.dbcfull;
	
	%outer
	udbcfullE = filtering.u.outer.dbcfull;
	vdbcfullE = filtering.v.outer.dbcfull;
	pdbcfullE = filtering.p.outer.dbcfull;
	qdbcfullE = filtering.q.outer.dbcfull;
	
	for i=1:numel(par.varnames)
		matlab.lang.makeValidName([par.varnames{i} 'dbcfull']);
		eval([par.varnames{i} 'dbcfull = filtering.' par.varnames{i} '.inner.dbcfull;']);
		eval([par.varnames{i} 'dbcfullE = filtering.' par.varnames{i} '.outer.dbcfull;']);
	end
	
	%gridify
	dbcfields = fields(udbcfull);
	
	for j=1:numel(dbcfields)
		
		Udbcfull.(dbcfields{j}) = gridify(udbcfull.(dbcfields{j}),nx-1,ny);
		Vdbcfull.(dbcfields{j}) = gridify(vdbcfull.(dbcfields{j}),nx,ny-1);
		Pdbcfull.(dbcfields{j}) = gridify(pdbcfull.(dbcfields{j}),nx,ny);
		Qdbcfull.(dbcfields{j}) = gridify(qdbcfull.(dbcfields{j}),nx-1,ny-1);
		
		UdbcfullE.(dbcfields{j}) = gridify(udbcfullE.(dbcfields{j}),nx+1,ny+2);
		VdbcfullE.(dbcfields{j}) = gridify(vdbcfullE.(dbcfields{j}),nx+2,ny+1);
		PdbcfullE.(dbcfields{j}) = gridify(pdbcfullE.(dbcfields{j}),nx+2,ny+2);
		QdbcfullE.(dbcfields{j}) = gridify(qdbcfullE.(dbcfields{j}),nx+1,ny+1);
		
	end
	
	uselEw = UdbcfullE.w&~(UdbcfullE.s|UdbcfullE.n);
	uselEe = UdbcfullE.e&~(UdbcfullE.s|UdbcfullE.n);
	%uselEs = UdbcfullE.s&~(UdbcfullE.w|UdbcfullE.e);
	%uselEn = UdbcfullE.n&~(UdbcfullE.w|UdbcfullE.e);	
	
	%vselEw = VdbcfullE.w&~(VdbcfullE.s|VdbcfullE.n);
	%vselEe = VdbcfullE.e&~(VdbcfullE.s|VdbcfullE.n);
	vselEs = VdbcfullE.s&~(VdbcfullE.w|VdbcfullE.e);
	vselEn = VdbcfullE.n&~(VdbcfullE.w|VdbcfullE.e);
	
	%inner boundary cond
	Urhs = reshape(filtering.u.inner.filterMat'*rhs.u.inner,[nx-1,ny])';
	Ubc = Urhs;
	Ubc(Udbcfull.w|Udbcfull.e) = par.dt/(hy^2*par.Re)*Urhs(Udbcfull.w|Udbcfull.e);
	Ubc(Udbcfull.n|Udbcfull.s) = 2*par.dt/(hx^2*par.Re)*Urhs(Udbcfull.n|Udbcfull.s);
	Ubc(Udbcfull.c) = 0;
	
	Vrhs = reshape(filtering.v.inner.filterMat'*rhs.v.inner*0,[nx,ny-1])';
	Vbc = Vrhs;
	Vbc(Vdbcfull.n|Vdbcfull.s) = par.dt/(hx^2*par.Re)*Vrhs(Vdbcfull.s|Vdbcfull.n);
	Vbc(Vdbcfull.w|Vdbcfull.e) = 2*par.dt/(hy^2*par.Re)*Vrhs(Vdbcfull.w|Vdbcfull.e);
	Vbc(Vdbcfull.c) = 0;
	
	QUbc = 0*reshape(filtering.q.inner.filterMat'*rhs.q.inner,[nx-1,ny-1])';
	QVbc = QUbc;
	
	%outer boundary cond
	UbcE = reshape(filtering.u.outer.filterMat'*rhs.u.outer,[nx+1,ny+2])';
	VbcE = reshape(filtering.v.outer.filterMat'*rhs.v.outer*0,[nx+2,ny+1])';
	QUbcE = 0*reshape(filtering.q.outer.filterMat'*rhs.q.outer,[nx+1,ny+1])';
	QVbcE = QUbcE;
	
	Lp = laplacian2(nx,ny,hx,hy,1,1,1,1,pdbcfull.w|pdbcfull.e,(pdbcfull.s|pdbcfull.n),filtering.p.inner.bciofull,-1);
	Lp = filtering.p.inner.filterMat*Lp*filtering.p.inner.filterMat';
	Lp(1,1) = 3/2*Lp(1,1);
	perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
	
	Lu = laplacian2(nx-1,ny,hx,hy,2,3,2,1,udbcfull.w|udbcfull.e,udbcfull.s|udbcfull.n,filtering.u.inner.bciofull,-1);
	Lu = par.dt/par.Re*Lu + speye(size(Lu,1));
	Lu = filtering.u.inner.filterMat*Lu*filtering.u.inner.filterMat';
	peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
	
	Lv = laplacian2(nx,ny-1,hx,hy,3,2,3,1,vdbcfull.w|vdbcfull.e,vdbcfull.s|vdbcfull.n,filtering.v.inner.bciofull,-1);
	Lv = par.dt/par.Re*Lv + speye(size(Lv,1));
	Lv = filtering.v.inner.filterMat*Lv*filtering.v.inner.filterMat';
	perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';
	
	Lq = laplacian2(nx-1,ny-1,hx,hy,2,2,2,1,qdbcfull.w|qdbcfull.e,qdbcfull.s|qdbcfull.n,filtering.q.inner.bciofull,-1);
	Lq = filtering.q.inner.filterMat*Lq*filtering.q.inner.filterMat';
	perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';
	
	if(par.useGPU)

		Rp = gpuArray(Lp);
		Rpt = Rp';

		Ru = gpuArray(Lu);
		Rut = Ru';

		Rv = gpuArray(Lv);
		Rvt = Rv';

		Rq = gpuArray(Lq);
		Rqt = Rq';

	end
	
	uValMat = reshape(filtering.u.inner.valind,[nx-1,ny])';
	uValMatE = reshape(filtering.u.outer.valind,[nx+1,ny+2])';
	uValMatInE = [zeros([ny,1]),uValMat,zeros([ny,1])];
	uValMatInE = logical([zeros([1,nx+1]);uValMatInE;zeros([1,nx+1])]);
	uvalindE = reshape(uValMatInE',[],1);
	
	vValMat = reshape(filtering.v.inner.valind,[nx,ny-1])';
	vValMatE = reshape(filtering.v.outer.valind,[nx+2,ny+1])';
	vValMatInE = [zeros([ny-1,1]),vValMat,zeros([ny-1,1])];
	vValMatInE = logical([zeros([1,nx+2]);vValMatInE;zeros([1,nx+2])]);
	vvalindE = reshape(vValMatInE',[],1);
	
	pValMat = reshape(filtering.p.inner.valind,[nx,ny])';
	pValMatE = reshape(filtering.p.outer.valind,[nx+2,ny+2])';
	pValMatInE = [zeros([ny,1]),pValMat,zeros([ny,1])];
	pValMatInE = logical([zeros([1,nx+2]);pValMatInE;zeros([1,nx+2])]);
	pvalindE = reshape(pValMatInE',[],1);
	
	qValMat = reshape(filtering.q.inner.valind,[nx-1,ny-1])';
	qValMatE = reshape(filtering.q.outer.valind,[nx+1,ny+1])';
	qValMatInE = [zeros([ny-1,1]),qValMat,zeros([ny-1,1])];
	qValMatInE = logical([zeros([1,nx+1]);qValMatInE;zeros([1,nx+1])]);
	qvalindE = reshape(qValMatInE',[],1);
	
	Pxsel = circshift(Pdbcfull.w,-1,2)&~circshift(reshape(filtering.p.inner.bciofull,[nx,ny])',-1,2);
	Pxsel = Pxsel(:,1:end-1);
	Pysel = Pdbcfull.n|circshift(Pdbcfull.s,-1,1);
	Pysel = Pysel(1:end-1,:);
	
	UselinN = circshift(UdbcfullE.n,-1)&~(grids.u.outer.Ymesh == min(grids.u.outer.yinit));
	UselinS = circshift(UdbcfullE.s,1)&~(grids.u.outer.Ymesh == max(grids.u.outer.yinit));
	VselinW = circshift(VdbcfullE.w,1,2)&~(grids.v.outer.Xmesh == max(grids.v.outer.xinit));
	VselinE = circshift(VdbcfullE.e,-1,2)&~(grids.v.outer.Xmesh == min(grids.v.outer.xinit));
	
	uvxsel = reshape(filtering.u.outer.valind,[nx+1,ny+2])';
	uvysel = reshape(filtering.v.outer.valind,[nx+2,ny+1])';
	
	fprintf('\n')
	lents = 0;
	for j=1:par.timesteps
		
		% nonlinear terms/explicit convection
		gamma = min(1.2*par.dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
		
		Ue(uValMatInE) = U(uValMat);
		Ue(uselEw|uselEe) = UbcE(uselEw|uselEe);
		Ue(UdbcfullE.n) = 2*UbcE(UdbcfullE.n)-Ue(UselinN);
		Ue(UdbcfullE.s) = 2*UbcE(UdbcfullE.s)-Ue(UselinS);
		
		Ve(vValMatInE) = V(vValMat);
		Ve(vselEs|vselEn) = VbcE(vselEs|vselEn);
		Ve(VdbcfullE.w) = 2*VbcE(VdbcfullE.w)-Ve(VselinW);
		Ve(VdbcfullE.e) = 2*VbcE(VdbcfullE.e)-Ve(VselinE);
		
		Ua = mvgavg(Ue);
		Ud = diff(Ue)/2;
		Va = mvgavg(Ve,2);
		Vd = diff(Ve')'/2;
		UVx = diff((Ua.*Va-gamma*abs(Ua).*Vd)')'/hx;
		UVy = diff((Ua.*Va-gamma*Ud.*abs(Va)))/hy;
		
		UVx(VdbcfullE.w(:,2:end-1)|VdbcfullE.e(:,2:end-1)|~uvysel(:,2:end-1)) = 0;
		UVy(UdbcfullE.s(2:end-1,:)|UdbcfullE.n(2:end-1,:)|~uvxsel(2:end-1,:)) = 0;
		
		Ue2 = Ue;
		Ue2(~uValMatInE&(UdbcfullE.n|UdbcfullE.s)) = 0;
		Ue2 = Ue2(2:end-1,:);
		Ve2 = Ve;
		Ve2(~vValMatInE&(VdbcfullE.w|VdbcfullE.e)) = 0;
		Ve2 = Ve2(:,2:end-1);
		Ua = mvgavg(Ue2,2);
		Ud = diff(Ue2')'/2;
		Va = mvgavg(Ve2);
		Vd = diff(Ve2)/2;
		U2x = diff((Ua.^2-gamma*abs(Ua).*Ud)')'/hx;
		V2y = diff((Va.^2-gamma*abs(Va).*Vd))/hy;
		
		U = U-par.dt*(UVy(:,2:end-1)+U2x);
		V = V-par.dt*(UVx(2:end-1,:)+V2y);
		
		% implicit diffusion
		rhsu = reshape((U+Ubc)',[],1);
		rhsu = rhsu(filtering.u.inner.valind);
		u(peru,1) = (Ru\(Rut\rhsu(peru)));
		U = reshape((1*filtering.u.inner.filterMat')*u,nx-1,ny)';
		
		rhsv = reshape((V+Vbc)',[],1);
		rhsv = rhsv(filtering.v.inner.valind);
		v(perv,1) = (Rv\(Rvt\rhsv(perv)));
		V = reshape((1*filtering.v.inner.filterMat')*v,nx,ny-1)';
		
		% pressure correction
		Uep = Ue;
		Uep(uValMatInE) = U(uValMat);
		Ue(uselEw|uselEe) = UbcE(uselEw|uselEe);
		Uep(UdbcfullE.s|UdbcfullE.n) = 0;
		Uep = Uep(2:end-1,:);
		
		Vep = Ve;
		Vep(vValMatInE) = V(vValMat);	
		Vep(vselEs|vselEn) = VbcE(vselEs|vselEn);
		Vep(VdbcfullE.w|VdbcfullE.e) = 0;
		Vep = Vep(:,2:end-1);
		
		rhsp = reshape((diff(Uep')'/hx+diff(Vep)/hy)',[],1);
		rhsp = rhsp(filtering.p.inner.valind);
		p(perp,1) = -(Rp\(Rpt\rhsp(perp)));
		P = reshape((1*filtering.p.inner.filterMat')*p,nx,ny)';
		Px = diff(P')'/hx;
		%Px(Pxsel) = 0;
		Py = diff(P)/hy;
		%Py(Pysel) = 0;
		U = U-Px;
		V = V-Py;
		
		Ue(uValMatInE) = U(uValMat);
		Ue(UdbcfullE.e) = UbcE(UdbcfullE.e);
		Ue(UdbcfullE.w) = UbcE(UdbcfullE.w);
		Ue(UdbcfullE.n) = UbcE(UdbcfullE.n);
		Ue(UdbcfullE.s) = UbcE(UdbcfullE.s);
		
		Ve(vValMatInE) = V(vValMat);
		Ve(VdbcfullE.s) = VbcE(VdbcfullE.s);
		Ve(VdbcfullE.n) = VbcE(VdbcfullE.n);
		Ve(VdbcfullE.w) = VbcE(VdbcfullE.w);
		Ve(VdbcfullE.e) = VbcE(VdbcfullE.e);
		
		if((j==1 || mod(j,par.plotoniter)==0)&&~par.noplot)
			
			%stream function
			%rhsq = reshape((diff(U)/hy-diff(V')'/hx)',[],1);
			%rhsq(filtering.q.inner.onfull) = rhs.q.inner(filtering.q.inner.on);
			%rhsq = rhsq(filtering.q.inner.valind);
			%q = Rq\(Rqt\rhsq(perq));
			%Q = reshape((1*filtering.q.inner.filterMat')*q,nx-1,ny-1)';
			
% 			res.U = gather(U);
% 			res.V = gather(V);
% 			res.P = gather(P);
% 
% 			res.Ue = gather(Ue);
% 			res.Ve = gather(Ve);

			[res.Ue,res.Ve] = makeQuiverData(Ue,Ve,QdbcfullE,QUbcE,QVbcE,reshape(filtering.q.outer.valind,[nx+1,ny+1])');
			res.U = zeros(ny-1,nx-1);
			res.V = zeros(ny-1,nx-1);
			res.U(qValMat) = res.Ue(qValMatInE);
			res.V(qValMat) = res.Ve(qValMatInE);
			
			a = res.Ue;
			b = res.U;
			res.P = gather(P);
			Qu = cumsum(res.U*hy);
			Qv = -cumsum(res.V'*hx)';
			res.Q = Qu;
			
			if(~exist('figs','var'))
				figs = PlotNS(grids,filtering,res,par);
			else
				figs = PlotNS(grids,filtering,res,par,figs);
			end
		end
		
		timestr = ['Time step: ' num2str(j*par.dt,'%5.2f') '/' num2str(par.tf,'%5.2f') 's'];
		fprintf([repmat('\b',1,lents) timestr]);
		lents = numel(timestr);
		
	end
	fprintf('\n');
end

function M = gridify(v,nx,ny)
	M = reshape(v,[nx,ny])';
end

function v = vectify(M)
	v = reshape(M',[],1);
end

function [Unew,Vnew] = makeQuiverData(U,V,QdbcfullE,QUbcE,QVbcE,qValMat)
	Unew = mvgavg(U);
	Vnew = mvgavg(V,2);
	Unew(QdbcfullE.s|QdbcfullE.n|QdbcfullE.c) = QUbcE(QdbcfullE.s|QdbcfullE.n|QdbcfullE.c);
	Vnew(QdbcfullE.w|QdbcfullE.e|QdbcfullE.c) = QVbcE(QdbcfullE.w|QdbcfullE.e|QdbcfullE.c);
	Unew(~(qValMat)) = 0;
	Vnew(~(qValMat)) = 0;
	%Len = sqrt(Unew.^2+Vnew.^2+eps);
	%Unew = Unew./Len;
	%Vnew = Vnew./Len;
end

