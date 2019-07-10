function histPlot = plotSampleHists(M,parameters,binCount) 
	
	for mo=1:M  

	    a = min(parameters(mo).sample);
	    b = max(parameters(mo).sample);  

	    bins=linspace(a,b,binCount);
	    
	    figure(1)

	    if M<=4
	        rm = ceil(M/2);
	        subplot(2,rm,mo)
	    elseif (4<M) && (M<=9)
	        rm = ceil(M/3);
	        subplot(3,rm,mo)
	    else
	        rm = ceil(M/4);
	        subplot(4,rm,mo)
	    end

	    histPlot = histogram(parameters(mo).sample,bins); 

	    hold on
	    box on

	    ylabel('# samples')

	    xlim([a b]);
	    xlabel('sampled values')

	    title(parameters(mo).name)  
	end %end mo loop through the parameters

end