function bool=effeq(a,b)
	%EFFEQ checks if numbers/arrays are effectively equal
	bool = abs(a-b) < eps;
end
