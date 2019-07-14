function prccPlot = plotVariedPRCC(M,N,x,labelstring,parameters,prcc)

    figure()
    hold on
    box on
    for mm=1:M
        plot(x, prcc(mm,:),'LineWidth',2.0);
    end
    
    prccPlot = gca; 

    xlabel('time (days)');
    ylabel('PRCC value');
    title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);
    legend(parameters.name,'Location','EastOutside')

end