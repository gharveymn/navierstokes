function figs = PlotNS(grids,filtering,res,par,figs)
	set(gcf, 'renderer', 'zbuffer')
	if(~exist('figs','var'))
		figs = InitialPlot(grids,filtering,res,par);
	else
		%update if we already initialized the plots
		figs = Update(grids,filtering,res,par,figs);
	end
end

function figs = InitialPlot(grids,filtering,res,par)
	
	U = res.U;
	V = res.V;
	P = res.P;
	Q = res.Q;

	figure(1)
	clf
	hold on
	ax = MakeAxis(grids.q.inner.Xmesh,grids.q.inner.Ymesh);
	%surf
	figs.f11 = surf(grids.q.inner.Xmesh,grids.q.inner.Ymesh,Q,'edgecolor','none','facecolor','interp');	
	%quiver
	[Unew,Vnew] = makeQuiverData(U,V);
	figs.f12 = quiver3(grids.q.inner.Xmesh,grids.q.inner.Ymesh,max(max(Q))*ones(size(Unew)),Unew,Vnew,zeros(size(Unew)));
	axis(ax)
	title('velocity vector field')
	hold off
	
	figure(2)
	clf
	ax = MakeAxis(grids.p.inner.Xmesh,grids.p.inner.Ymesh);
	figs.f2 = surf(grids.p.inner.Xmesh,grids.p.inner.Ymesh,P);
	axis(ax)
	title('$p$','interpreter','latex','FontSize',20)

	%contour

	figure(3)
	clf
	[C1,h1] = contour(grids.u.inner.Xmesh,grids.u.inner.Ymesh,U,par.conlines);
	figs.f3 = {C1,h1};
	axis(ax)
	title('$u$','interpreter','latex','FontSize',20)

	figure(4)
	clf
	[C2,h2] = contour(grids.v.inner.Xmesh,grids.v.inner.Ymesh,V,par.conlines);
	figs.f4 = {C2,h2};
	axis(ax)
	title('$v$','interpreter','latex','FontSize',20)

	figure(5)
	clf
	[C3,h3] = contour(grids.q.inner.Xmesh,grids.q.inner.Ymesh,Q,par.conlines);
	figs.f5 = {C3,h3};
	axis(ax)
	title('$\psi$','interpreter','latex','FontSize',20)
	
	drawnow;
	
end

function figs = Update(grids,filtering,res,par,figs)
	
	lastwarn('')
	
	try
		U = res.U;
		V = res.V;
		P = res.P;
		Q = res.Q;
		[Unew,Vnew] = makeQuiverData(U,V);

		set(figs.f11,'ZData',Q);
		set(figs.f12,'ZData',max(max(Q))*ones(size(Unew)));
		set(figs.f12,'UData',Unew);
		set(figs.f12,'VData',Vnew);
		set(figs.f2,'ZData',P);
		set(figs.f3{2},'ZData',U);
		set(figs.f4{2},'ZData',V);
		set(figs.f5{2},'ZData',Q);

		drawnow;
		
	catch ME
		disp('Couldn''t update one of the figures; we''ll try to make new ones')
		figs = InitialPlot(grids,filtering,res,par);
	end
	
	if(~isempty(lastwarn))
		figs = InitialPlot(grids,filtering,res,par);
	end
	
end

function ax = MakeAxis(x,y)
	x = x(isfinite(x));
	y = y(isfinite(y));
	
	minx = min(min(x));
	maxx = max(max(x));
	miny = min(min(y));
	maxy = max(max(y));
	absmin = min(minx,miny);
	absmax = max(maxx,maxy);
	centerx = (maxx+minx)/2;
	centery = (maxy+miny)/2;
	difa = absmax-absmin;
	ax = [centerx-difa/2, centerx+difa/2, centery-difa/2, centery+difa/2];
end

function clrs = MakeClrs(v)
	clr = abs(v)./norm(v(isfinite(v)),inf);
	clrs = [clr zeros(numel(clr),1) 1-clr];
end

function [Unew,Vnew] = makeQuiverData(U,V)
	Unew = mvgavg(U);
	Vnew = mvgavg(V,2);
	Len = sqrt(Unew.^2+Vnew.^2+eps);
	Unew = Unew./Len;
	Vnew = Vnew./Len;
end