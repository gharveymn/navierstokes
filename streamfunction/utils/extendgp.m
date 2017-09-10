function rhs = extendgp(rhs,bcfull,valindouter,gpca,nx)
	%EXTENDGP extends rhs to ghostpoints order
	%   Detailed explanation goes here
	
	bcw = bcfull{1};
	bce = bcfull{2};
	bcs = bcfull{3};
	bcn = bcfull{4};
	bcc = bcfull{5};
	
	%set outer regions
	
	%corners (up to second order)
	rhs = makecorners(rhs,bcc,valindouter,gpca,nx);
	numIter = numel(gpca);
	
	for i=1:numIter
		%west
		bcwi = circshift(bcw,-i);
		rhs(bcwi(valindouter)) = rhs(bcw(valindouter));
		
		%east
		bcei = circshift(bce,i);
		rhs(bcei(valindouter)) = rhs(bce(valindouter));
		
		%south
		bcsi = circshift(bcs,-i*nx);
		rhs(bcsi(valindouter)) = rhs(bcs(valindouter));
		
		%north
		bcni = circshift(bcn,i*nx);
		rhs(bcni(valindouter)) = rhs(bcn(valindouter));
	end
	
end

function rhs = makecorners(rhs,bcc,valind,gpca,nx)
	%input i for the order
	
	for k=1:numel(gpca)
		gp = gpca{k};
		for i=-k:k
			for j=-k:k
				bccsw = circshift(circshift(bcc,-i),-j*nx) & gp;
				bccswr = circshift(circshift(bccsw,j*nx),i);
				
				bccnw = circshift(circshift(bcc,-i),j*nx) & gp;
				bccnwr = circshift(circshift(bccnw,-j*nx),i);
				
				bccse = circshift(circshift(bcc,i),-j*nx) & gp;
				bccser = circshift(circshift(bccse,j*nx),-i);
				
				bccne = circshift(circshift(bcc,i),j*nx) & gp;
				bccner = circshift(circshift(bccne,-j*nx),-i);
				
				rhs(bccsw(valind)) = rhs(bccswr(valind));
				rhs(bccnw(valind)) = rhs(bccnwr(valind));
				rhs(bccse(valind)) = rhs(bccser(valind));
				rhs(bccne(valind)) = rhs(bccner(valind));
			end
		end
	end
end

