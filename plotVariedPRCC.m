function prccPlot = plotVariedPRCC(M,N,x,labelstring,parameters,prcc,studentT)

    figure()
    hold on
    box on
    for mm=1:M
        plot(x, prcc(mm,:))
    end
    
    prccPlot = gca; 

    xlabel('time (days)');
    ylabel('PRCC value');
    title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);
    legend(parameters.name,'Location','EastOutside')



    if N>250 % plot StudentT
        figure()
        hold on
        box on

        SigPt05=linspace(1.658,1.658,length(x)); %These levels of significance are true for N = 120 or more
        SigPt01=linspace(2.358,2.358,length(x));
        SigPt001=linspace(3.373,3.373,length(x));

        for mn=1:M
            plot(x, studentT(mn,:))
        end

        plot(x, SigPt05, '-.b');
        plot(x, SigPt01, '-.r');
        plot(x, SigPt001, '-.k');
        plot(x, -SigPt05, '-.b');
        plot(x, -SigPt01, '-.r');
        plot(x, -SigPt001, '-.k');

        xlabel('time (days)');
        ylabel('t value');% = PRCC*sqrt((Number of Samples-2)/(1-PRCC))');
        legend(parameters.name,'Location','EastOutside')
        title(['Sensitivity Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples, at the center']);

        figurelabel1=(['LHSgeneral-N',num2str(N),'-', labelstring,'-StudentT.fig']); %assign an appropriate filename
        figurelabel2=(['LHSgeneral-N',num2str(N),'-', labelstring,'-StudentT.png']);
        pause(5) %Time to dock/maximize the figure before it saves, if you prefer.
        saveas(gcf,figurelabel1); %save fig
        saveas(gcf,figurelabel2); %save png
    end;
end