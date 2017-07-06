function figs = Plot(mat,vec,par,figs)
	
	if(nargin == 3)
		figs = InitialPlot(mat,vec,par);
	else
		%update if we already initialized the plots
		figs = Update(mat,vec,par,figs);
	end
end

function figs = InitialPlot(mat,vec,par)
	
	X = mat(:,:,1);
	Y = mat(:,:,2);
	U = mat(:,:,3);
	V = mat(:,:,4);
	Q = mat(:,:,5);
	
	x = vec(:,1);
	y = vec(:,2);
	u = vec(:,3);
	v = vec(:,4);
	q = vec(:,5);
	
	ax = MakeAxis(x,y);
	
	if(par.toPlot == 1)
		%surf
		
		figure(1)
		f1 = surf(X,Y,U,'edgecolor','none','facecolor','interp');
		axis(ax)
		title('$u$','interpreter','latex','FontSize',20)

		figure(2)
		f2 = surf(X,Y,V,'edgecolor','none','facecolor','interp');
		axis(ax)
		title('$v$','interpreter','latex','FontSize',20)
		
		figure(3)
		f3 = surf(X,Y,Q,'edgecolor','none','facecolor','interp');
		%f3 = contour(X,Y,Q);
		axis(ax)
		title('$\psi$','interpreter','latex','FontSize',20)
		
		figs = {f1,f2,f3};
		
	elseif(par.toPlot == 2)
		%quiver
		
		figure(1)
		ax = MakeAxis(X,Y);
		f1 = quiver(x,y,u,v);
		axis(ax)
		title('velocity vector field')

		figure(2)
		ax = MakeAxis(X,Y);
		clrs = abs(q)./max(abs(q(isfinite(q))));
		f2 = scatter3(x,y,q,10,[clrs zeros(numel(clrs),1) 1-clrs], '.');
		axis(ax)
		title('$\psi$','interpreter','latex','FontSize',20)
		
		figs = {f1,f2};
		
	elseif(par.toPlot == 3)
		%scatter
		
		figure(1)
		clrs = MakeClrs(u);
		f1 = scatter3(x,y,u,10,clrs,'.');
		axis(ax)
		title('$u$','interpreter','latex','FontSize',20)

		figure(2)
		clrs = MakeClrs(v);
		f2 = scatter3(x,y,v,10,clrs,'.');
		axis(ax)
		title('$v$','interpreter','latex','FontSize',20)

		figure(3)
		clrs = MakeClrs(q);
		f3 = scatter3(x,y,q,10,clrs,'.');
		axis(ax)
		title('$\psi$','interpreter','latex','FontSize',20)
		
		figs = {f1,f2,f3};
		
	elseif(par.toPlot == 4)
		%contour
		
		figure(1)
		[C1,h1] = contour(X,Y,U,par.conlines);
		f1 = {C1,h1};
		axis(ax)
		title('$u$','interpreter','latex','FontSize',20)

		figure(2)
		[C2,h2] = contour(X,Y,V,par.conlines);
		f2 = {C2,h2};
		axis(ax)
		title('$v$','interpreter','latex','FontSize',20)
		
		figure(3)
		[C3,h3] = contour(X,Y,Q,par.conlines);
		f3 = {C3,h3};
		axis(ax)
		title('$\psi$','interpreter','latex','FontSize',20)
		
		figs = {f1,f2,f3};
		
	end
	
	drawnow;
	
end

function figs = Update(mat,vec,par,figs)
	
	lastwarn('')
	
	try
		U = mat(:,:,3);
		V = mat(:,:,4);
		Q = mat(:,:,5);

		u = vec(:,3);
		v = vec(:,4);
		q = vec(:,5);

		if(par.toPlot == 1)

			set(figs{1},'ZData',U);
			set(figs{2},'ZData',V);
			set(figs{3},'ZData',Q);

		elseif(par.toPlot == 2)

			set(figs{1},'UData',u);
			set(figs{1},'VData',v);

			clrs = MakeClrs(q);
			set(figs{2},'ZData',q);
			set(figs{2},'CData',clrs);

		elseif(par.toPlot == 3)

			clrs = MakeClrs(u);
			set(figs{1},'ZData',u);
			set(figs{1},'CData',clrs);

			clrs = MakeClrs(v);
			set(figs{2},'ZData',v);
			set(figs{2},'CData',clrs);

			clrs = MakeClrs(q);
			set(figs{3},'ZData',q);
			set(figs{3},'CData',clrs);
			
		elseif(par.toPlot == 4)
			set(figs{1}{2},'ZData',U);
			set(figs{2}{2},'ZData',V);
			set(figs{3}{2},'ZData',Q);
		end

		drawnow;
		
	catch ME
		disp('Couldn''t update one of the figures; we''ll try to make new ones')
		figs = InitialPlot(mat,vec,par);
	end
	
	if(~isempty(lastwarn))
		figs = InitialPlot(mat,vec,par);
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



