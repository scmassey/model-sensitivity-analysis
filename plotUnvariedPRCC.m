function prccPlot = plotUnvariedPRCC(M,N,x,labelstring,parameters,prcc,studentT)

    figure()
    hold on
    box on
    
    p=categorical({parameters.name});
    prccPlot = bar(p,prcc(1:M));

    % xlabel('sampled parameters');
    ylabel('PRCC value');

    title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);

end