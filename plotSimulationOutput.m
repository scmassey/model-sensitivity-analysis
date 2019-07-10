function outPlot = plotSimulationOutput(x,N,Simdata,labelstring) 

	TL=length(x);
	%Creating L, matrix of all outputs at each time, for each of the N samples
	L=zeros(TL,N);
	R=zeros(TL,N);  
	for i=1:TL
	    for k=1:N
	            L(i,k)=Simdata(k).y(i);
	    end
	end
	Z=L';% {want mean and std dev computed with each day kept together, looking at variability across samples,
	     % {and returned as a row vector rather than a column vector (why use Z=L' instead of L)
	pause(.2)

	figure() 

	Y=mean(Z);   
	E=std(Z);

	outPlot = errorbar(x,Y,E)

	xlabel('x')

	ylabel(labelstring)

	title(['Error bar plot of ',labelstring,' from LHS simulations']);

end
	