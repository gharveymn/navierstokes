function figs = Plot(grids,filtering,res,par,figs)
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
	Q = res.Q;
	ax = MakeAxis(grids.inner.Xmesh,grids.inner.Ymesh);

	if(par.toPlot == 1)
		
		figure(1)
		clf
		hold on
		%surf
		figs.f11 = surf(grids.outer.Xmesh,grids.outer.Ymesh,Q,'edgecolor','none','facecolor','interp');
		%quiver
		figs.f12 = quiver3(grids.outer.Xmesh,grids.outer.Ymesh,max(max(Q))*ones(size(U)),U,V,zeros(size(U)),par.quivVectSca,'k-');
		axis(ax)
		title('velocity vector field')
		hold off

		%contour

		figure(2)
		clf
		[C1,h1] = contour(grids.outer.Xmesh,grids.outer.Ymesh,U,par.conlines);
		figs.f2 = {C1,h1};
		axis(ax)
		title('$u$','interpreter','latex','FontSize',20)

		figure(3)
		clf
		[C2,h2] = contour(grids.outer.Xmesh,grids.outer.Ymesh,V,par.conlines);
		figs.f3 = {C2,h2};
		axis(ax)
		title('$v$','interpreter','latex','FontSize',20)

		figure(4)
		clf
		[C3,h3] = contour(grids.outer.Xmesh,grids.outer.Ymesh,Q,par.conlines);
		figs.f4 = {C3,h3};
		axis(ax)
		title('$\psi$','interpreter','latex','FontSize',20)
		
	elseif(par.toPlot == 2)
		
		figure(1)
		clf
		figs.f1 = surf(grids.outer.Xmesh,grids.outer.Ymesh,U);
		axis(ax)
		title('$U$','interpreter','latex','FontSize',20)
		
		figure(2)
		clf
		figs.f2 = surf(grids.outer.Xmesh,grids.outer.Ymesh,V);
		axis(ax)
		title('$V$','interpreter','latex','FontSize',20)
		
		figure(3)
		clf
		figs.f3 = surf(grids.outer.Xmesh,grids.outer.Ymesh,Q);
		axis(ax)
		title('$Q$','interpreter','latex','FontSize',20)
	
	end
	
	drawnow;
	
	figure(1)
	
end

function figs = Update(grids,filtering,res,par,figs)
	
	lastwarn('')
	
	try
		U = res.U;
		V = res.V;
		Q = res.Q;
		
		if(par.toPlot == 1)

			set(figs.f11,'ZData',Q);
			set(figs.f12,'ZData',max(max(Q))*ones(size(U)));
			set(figs.f12,'UData',U);
			set(figs.f12,'VData',V);
			umax = max(max(U));
			umin = min(min(U));
			
			vmax = max(max(V));
			vmin = min(min(V));
			
% 			pmax = max(max(P));
% 			pmin = min(min(P));
			
			qmax = max(max(Q));
			qmin = min(min(Q));
			
			set(figs.f2{2},'LevelList',linspace(umax,umin,par.conlines));
 			set(figs.f2{2},'ZData',U);
			
			set(figs.f3{2},'LevelList',linspace(vmax,vmin,par.conlines));
 			set(figs.f3{2},'ZData',V);
			
			set(figs.f4{2},'LevelList',linspace(qmax,qmin,par.conlines));
 			set(figs.f4{2},'ZData',Q);

% 			ax = MakeAxis(grids.q.inner.Xmesh,grids.q.inner.Ymesh);
% 
% 			figure(3)
% 			[C1,h1] = contour(grids.u.inner.Xmesh,grids.u.inner.Ymesh,U,par.conlines);
% 			figs.f3 = {C1,h1};
% 			axis(ax)
% 			title('$u$','interpreter','latex','FontSize',20)
% 
% 			figure(4)
% 			[C2,h2] = contour(grids.v.inner.Xmesh,grids.v.inner.Ymesh,V,par.conlines);
% 			figs.f4 = {C2,h2};
% 			axis(ax)
% 			title('$v$','interpreter','latex','FontSize',20)
% 
% 			figure(5)
% 			[C3,h3] = contour(grids.q.inner.Xmesh,grids.q.inner.Ymesh,Q,par.conlines);
% 			figs.f5 = {C3,h3};
% 			axis(ax)
% 			title('$\psi$','interpreter','latex','FontSize',20)

		elseif(par.toPlot == 2)
			set(figs.f1,'ZData',U);
			set(figs.f2,'ZData',V);
			set(figs.f3,'ZData',Q);
		end

		drawnow;
		
	catch ME
		%throw(ME)
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