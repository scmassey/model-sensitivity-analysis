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

% TO DO: Add save command here so that as alternative, 
% can prompt to load previously sampled parameters 
% (in case want to re-run with diff time limits?)

%---Plot Histograms to Visualize Parameter Distributions---% 
% TO DO: ask for feedback - nest this into the drawsamples function?                                                                                      % 

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

% EDIT THE FOLLOWING FILE NAME FOR SAVING YOUR SIMULATION DATA:
save LHS_testlinear   % give the workspace an appropriate and unique name


% TO DO: ask for feedback - make this a separate function call? 
%---Visualize the spread of simulations results with errorbar plots---%
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


%__________________________________________________________________________
%
% THIRD STEP: Compute the Partial Rank Correlation Coefficients (PRCC)
%__________________________________________________________________________

K=zeros(N,M+1); % number of samples drawn N x number of model parameters (M) + 1 output column
K([1:N],[1:M])=A([1:N],[1:M]); %Taking the information from A and putting it into K

% EDIT THE FOLLOWING STRING TO NAME YOUR SIMULATION DATA:
labelstring = 'y' 
    
tx_varied = input('Do you want to look at results that might change over time or space (''Y'' or ''N'')?: ');

if strcmp(tx_varied,'Y')

% EDIT INDICES FOR EVALUATING RESULTS AS NEEDED:
    % To look at variation in TIME: 
    % x_index = 5; % specify spatial index if a PDE model
    % x_index = []; % for time varied ODE model
%     [prcc,studentT]=TimeVariedPRCC(M,N,K,Simdata,t,x_index);

    % To look at variation in SPACE:
    % t_index= 10; % specify time index if a PDE model
    t_index =[]; % for a spatially varied ODE model
    [prcc,studentT]=VariedPRCC(M,N,K,Simdata,x,t_index);
    
elseif strcmp(tx_varied,'N')

    % EDIT THE FOLLOWING INDICES FOR EVALUATING RESULTS:
    [x_index,t_index]=[1,]; %if an ODE, leave the nonapplicable index empty; if PDE, specify indices for both
    [prcc,studentT]=UnvariedPRCC(M,N,K,Simdata,t_index,x_index);

end

% TO DO: change plot code to do waterfalls for unvaried - and move to that function call?

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
