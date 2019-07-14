function outPlot = plotSimulationOutput(x,N,OutputOfInterest,labelstring) 

	figure() 

	Y=mean(OutputOfInterest);   
	E=std(OutputOfInterest);

	outPlot = errorbar(x,Y,E);

	xlabel('x')

	ylabel(labelstring)

	title(['Error bar plot of ',labelstring,' from LHS simulations']);

end
	