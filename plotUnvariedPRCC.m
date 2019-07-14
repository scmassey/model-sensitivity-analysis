function prccPlot = plotUnvariedPRCC(M,N,x,labelstring,parameters,prcc)

    figure()
    hold on
    box on

    % parNames = zeros(1,N);

    % for i=1:N
    % 	parNames(i) = getfield(parameters(i),'name');
    % end 
    % p=categorical({parNames})

    temp = {parameters.name}
    p=categorical(temp);

    prccPlot = bar(p,prcc(1:M));

    % xlabel('sampled parameters');
    ylabel('PRCC value');

    title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);

end