function [x2,R] = maxnormalize(x)
	x2 = (x-min(x))./(max(x)-min(x));
	R = range(x);
	
