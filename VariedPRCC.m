function [prcc,studentT]=VariedPRCC(M,N,K,Simdata,t,idx)

R=zeros(N,M+1); % initializing the rank matrix

prcc=zeros(M+1,length(t)); % initializing matrix to hold PRCC statistic for each parameter at each time

if N>=250
    studentT=zeros(M+1,length(t)); % initializing matrix to hold the studentT statistic for each parameter at each time
else
    studentT = [];
    fprintf(['Not enough samples (<250) for computing studentT statistic. \n'])
end;


for steps=1:length(t)  % loop through time (or space) - TO DO: would like to be able to do both and get surfaces..

    for n=1:N %loop through each of the N samples
%             output=Simdata(n).r(steps,space)/(Simdata(n).r(steps,space)+Simdata(n).c(steps,space));%output variable of interest (can be a ratio of outputs)
        output=Simdata(n).y(steps);
        K(n,M+1)=output; %Output results go in the last column, the M+1st column.
    end

    for m=1:M+1 %for each parameter rank the data, including the output data
        for n=1:N
%          select parameter column m
        [s,i]=sort(K(:,m));  %sort according to the rank of column m 
        r(i,1)=[1:N]'; 
        R(n,m)=r(n,1); %store the ranking of the parameter at it's position in K
        end
    end

    C=corrcoef(R);
    if det(C)<=10^-16 %If the determinant is singular or very nearly singular
        B=pinv(C);    %must use pseudo inverse (as inv will not be accurate).
        ST=num2str(steps); %Report the timestep of the singularity (since space is fixed).
        fprintf(['C is singular at timestep ', ST,'. \n'])
    else
        B=inv(C); %The dterminant is not singular, so the inverse is valid.
    end

    if N>=250
        for w=1:M %iterate thru the M parameters to calculate PRCCs for each with ratio
            prcc(w,steps)=(-B(w,M+1))/sqrt(B(w,w)*B(M+1,M+1)); % the PRCC between the wth parameter and the ratio
            studentT(w,steps)=prcc(w,steps)*sqrt((N-2)/(1-prcc(w,steps))); % the studentT statistic corresponding to the PRCC
        end
    else
        for w=1:M %iterate thru the M parameters to calculate PRCCs for each with ratio
            prcc(w,steps)=(-B(w,M+1))/sqrt(B(w,w)*B(M+1,M+1)); % the PRCC between the wth parameter and the ratio
        end
    end;
    %Each studentT is the distribution (showing the significance of) the corresponding gamma/PRCC. 
end; 

%---Plot PRCCs---%
figure(1)
hold on
box on
for mm=1:M
    plot(x, prcc(mm,:))
end

xlabel('time (days)');
ylabel('PRCC value');
title(['Correlation Plot of ',labelstring,' from LHS simulations, ',num2str(N),' samples']);
legend(parameters.name,'Location','EastOutside')

figurelabel1=(['LHS-N',num2str(N),'-', labelstring,'-PRCC.fig']); %assign an appropriate filename
figurelabel2=(['LHS-N',num2str(N),'-', labelstring,'-PRCC.png']);
pause(5) %Time to dock/maximize the figure before it saves, if you prefer.
saveas(gcf,figurelabel1); %save fig
saveas(gcf,figurelabel2); %save png

if N>250 % plot StudentT
    figure(2)
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
