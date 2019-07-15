function prccPlot = plotUnvariedPRCC(M,N,labelstring,parameters,prcc)

    figure()
    hold on
    box on

    p=categorical({parameters.name});

    prccPlot = bar(p,prcc(1:M));

    ylabel('PRCC value');

    title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);

end