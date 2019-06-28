%Susan Massey
%Most General version of the LHS code
%Modified Sept 2018 and then Jun 2019 to be more interactive/usable

clear all
close all
clc

%__________________________________________________________________________
%
% FIRST STEP: Latin Hypercube Sampling of Model Parameters 
%__________________________________________________________________________

%User Input Needed Here:
M = input('Number of parameters to sample?: '); % Don't include add'l parameters you would like to leave fixed
while rem(M,1)~=0 || M<=0
    M = input('Number of parameters should be an integer > 0. Please re-enter the number of parameters: ');
end

N = input('Number of samples to draw? (recommend 25 to 1000; Higher numbers yield better results, but also takes longer): '); 
while rem(N,1)~=0 || N<=0
    N = input('Number of parameters should be an integer > 0. Please re-enter the number of parameters: ');
end

[parameters,A] = drawsamples(M,N);

%__________________________________________________________________________
%
% SECOND STEP: Solve Model Function with Sampled Parameters
%__________________________________________________________________________

% EDIT THE FOLLOWING VARIABLES, UNSAMPLED PARAMETERS, AND ANY OTHER ARGS HERE

% Specify the independent variables you will pass to the model function in a vector
x = linspace(0,10,101); % spatial domain for my testlinear function
% t = linspace(0,17,171);% time vector w/ number of time steps desired for kineticpdgfpdes

% Specify values of any unsampled parameters to pass to the function
unsampledps =  2; 

% Don't forget any solver parameters if needed: 
% w = 2; %spherical symmetry = 2 for pdepe.m to solve kineticpdgfpdes

Simdata = struct; %place where all output data will be stored for computing PRCCs 

tic % start measuring time to solve equations
for j=1:N   
    sampledparams=A(j,1:M);
    fprintf('Parameters passed for sample: ');fprintf('%u',j);fprintf(' of ');fprintf('%u\n',N);
    
%  EDIT THE FOLLOWING FUNCTION CALL--can name multiple outputs as in commented example below
%  Also make sure everything passes to your function correctly and all initial conditions,
%  etc. are specified before entering this for loop over j.
    Simdata(j).y=testlinear(x,sampledparams,unsampledps);
% %     sol = pdepe(w,@kineticpdgfpdes,@pdgfpdesic,@pdgfpdesbc,x,t); %make sure the right set of pdes are used
% %     Simdata(j).c = sol(:,:,1); % assign these appropriately to store the data from your solution vector in the struct
% %     Simdata(j).p = sol(:,:,2); % For example, the part of the solution corresponding to the c pde I called Simdata(j).c
% %     Simdata(j).r = sol(:,:,3); % and the (j) index is required to hold each sample. To call a specific point:
% %     Simdata(j).n = sol(:,:,4); % Simdata(j=3).c(timestep=25,spacestep=50), or whatever makes sense for your equations.
    toc
end

% EDIT THE FOLLOWING FILE NAME FOR SAVING YOUR SIMULATION DATA.
save LHS_testlinear %give the workspace an appropriate unique name

%__________________________________________________________________________
%
% THIRD STEP: Compute the Partial Rank Correlation Coefficients (PRCC)
%__________________________________________________________________________

K=zeros(N,M+1); % number of samples drawn N x number of model parameters (M) + 1 output column
K([1:N],[1:M])=A([1:N],[1:M]); %Taking the information from A and putting it into K

% Name for the model output in plots and filenames?
labelstring = 'y' 
    
time_varied = input('Do you want to look at results that might change over time or space (''Y'' or ''N'')?: ');

if strcmp(time_varied,'Y')

    % if also spatially varied, specify location index 
    % loc = 5;
    % vice versa if looking at spatial variation and need to fix timepoint:
    % time = 10;

    [prcc,studentT]=TimeVariedPRCC(M,N,K,Simdata,t,loc)
    % use [prcc,studentT]=TimeVariedPRCC(M,N,K,Simdata,x,time) if opposite
    
elseif strcmp(time_varied,'N')
    x=1;
    % t=5;
    % [x,t]=[1,5];
    
    prcc=zeros(M+1,1); % will hold PRCC statistic for each parameter for chosen point (x=0, or otherwise)
    studentT=zeros(M+1,1); % will hold the studentT statistic
    
    for n=1:N %loop through each of the N samples
        output=Simdata(n).y(x);
        K(n,M+1)=output; %Output results go in the added last column, the M+1st column.
    end

    for ml=1:M+1 %for each parameter rank the data, including the output data
        for nl=1:N
%          select parameter column m
        [s,i]=sort(K(:,ml));  %sort according to the rank of column m 
        r(i,1)=[1:N]'; 
        R(nl,ml)=r(nl,1); %store the ranking of the parameter at it's position in K
        end
    end

    C=corrcoef(R);
    if det(C)<=10^-16 %If the determinant is singular or very nearly singular
        B=pinv(C);    %must use pseudo inverse (as inv will not be accurate).
        fprintf(['C is singular, used pseudoinverse. \n'])
    else
        B=inv(C); %The determinant is not singular, so the inverse is valid.
    end

    for w=1:M %iterate thru the M parameters to calculate PRCCs for each with ratio
        prcc(w)=(-B(w,M+1))/sqrt(B(w,w)*B(M+1,M+1)); % the PRCC between the wth parameter and the ratio
        studentT(w)=prcc(w)*sqrt((N-2)/(1-prcc(w))); % the studentT statistic corresponding to the PRCC
    end
    %Each studentT is the distribution (showing the significance of) the
    %corresponding gamma/PRCC
end

%__________________________________________________________________________
%
% FOURTH STEP: Generate Plots to Display the PRCC Results
%__________________________________________________________________________

% EDIT THE FOLLOWING - IN EACH PLOT, MAKE SURE THE INDEPENDENT VARIABLE IS
% SPECIFIED CORRECTLY. (The variable x in each should maybe be renamed
% depending on what it was called above before running model simulations.)

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


%---Plot StudentTs---%
figure(2)
hold on
box on

%TO DO: add a check on N and only display if enough samples - if not, add
%an fprintf to tell the user..
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

%---Plot Histograms---%                                                                                       

% % mid = [mu(1),       mu(2),       mu(3),       mu(4),       mu(5),       mu(6),       mu(7),       mu(8),       mu(9),       mean(K(:,M))];
% % low = [a_vector(1), a_vector(2), a_vector(3), a_vector(4), a_vector(5), a_vector(6), a_vector(7), a_vector(8), a_vector(9), min(K(:,M))];  
% % high= [b_vector(1), b_vector(2), b_vector(3), b_vector(4), b_vector(5), b_vector(6), b_vector(7), b_vector(8), b_vector(9), max(K(:,M))];

for mo=1:M  
    a = min(parameters(mo).sample);
    b = max(parameters(mo).sample);                                                                      
    bins=linspace(a,b,10);   %for now, added the calculated highs and lows of the nondim params to the vectors 
    figure(5)
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
    % % % TO DO: make the subplot call dynamic
    % % subplot(4,3,m) % Make this the size needed for the number M of parameters you have sampled.
    hist(K(:,mo),bins);
    hold on
    box on
    ylabel('# samples')
    xlim([a b]);
    xlabel('sampled values')
    title(parameters(mo).name)  
end %end mo loop through the parameter columns, including nondimensional ones
pause(5) %Time to dock/maximize the figure before it saves, if you prefer.
figurelabel1=(['LHS-N',num2str(N),'-histograms.fig']);
figurelabel2=(['LHS-N',num2str(N),'-histograms.png']);
saveas(gcf, figurelabel1);
saveas(gcf, figurelabel2);

%---Get Data for Errorbar and Box Plots---%
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

%---Plot Errorbars---%
figure(6)       
Y=mean(Z);   
E=std(Z);
errorbar(x,Y,E)
xlabel('x')
ylabel(labelstring)
title(['Error bar plot of ',labelstring,' from LHS simulations']);
pause(5)
figurelabel1=(['LHSgeneral-N',num2str(N),'-',labelstring,'-ErrorbarPlot.fig']);
figurelabel2=(['LHSgeneral-N',num2str(N),'-',labelstring,'-ErrorbarPlot.png']);
saveas(gcf,figurelabel1);
saveas(gcf,figurelabel2);

% % % %---Plot Boxplots---%  This NEEDS BOX PLOT CODES, found in the general LHS code folder (the ones that Gargi found online)
% % % ZZ=zeros(N,length(x));
% % % for i=1:length(x)
% % %     for j=1:N
% % %         ZZ(j,i)=Z(j,i*10);
% % %     end
% % % end
% % % 
% % % figure(7)
% % % boxplotC(ZZ,0,'*',1,1.5,'b',0,1);
% % % xlabel('x')
% % % ylabel(labelstring)
% % % title(['Box plot of ',labelstring,' from LHSgeneral simulations at the center']);
% % % pause(5)
% % % figurelabel1=(['LHSgeneral-N',num2str(N),'-', labelstring,'-BoxPlot.fig']);
% % % figurelabel2=(['LHSgeneral-N',num2str(N),'-', labelstring,'-BoxPlot.png']);
% % % saveas(gcf,figurelabel1);
% % % saveas(gcf,figurelabel2);