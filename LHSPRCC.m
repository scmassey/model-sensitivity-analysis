%Susan Massey
%Most General version of the LHS code
%Modified Sept 2018 and then Jun 2019 to be more interactive/usable

clear all
close all
clc

%%_________________________________________________________________________
%
% FIRST STEP: Latin Hypercube Sampling of Model Parameters 
%__________________________________________________________________________

% This requests user input via the command line to get started:
M = input('Number of parameters to sample?: '); % Don't include add'l parameters you would like to leave fixed
while rem(M,1)~=0 || M<=0
    M = input('Number of parameters should be an integer > 0. Please re-enter the number of parameters: ');
end

N = input('Number of samples to draw? (recommend 100 to 1000): '); 
while rem(N,1)~=0 || N<=0
    N = input('Number of samples should be an integer > 0. Please re-enter the number of samples: ');
end

% This code will prompt for sample distribution specifics and return drawn and randomly paired parameter samples
[parameters,A] = DrawSamples(M,N);

% EDIT THE FOLLOWING FOR SAVING YOUR SAMPLES AND SIMULATION DATA:
outFileStr = 'LHS-testlinear'; % give workspace an appropriate unique name
outFileName1 = [outFileStr,'_samples.mat'];
outFileName2 = [outFileStr,'_results.mat'];

save(outFileName1, 'parameters', 'A')


%---Plot Histograms to Visualize Parameter Distributions---% 
 
% Specify how many bins you want for your histograms:                            
binCount = ceil(N/10); % can change this manually, just make sure it's an integer value  

histPlot = plotSampleHists(M,parameters,binCount);

pause(1) %Time to dock/maximize the figure before it saves, if you prefer.

% Name and save histograms visualizing the sample distributions:
figurelabel1=([outFileStr,'-N',num2str(N),'-histograms.fig']);
figurelabel2=([outFileStr,'-N',num2str(N),'-histograms.png']);
saveas(histPlot, figurelabel1);
saveas(histPlot, figurelabel2);

%%_________________________________________________________________________
%
% SECOND STEP: Solve Model Function with Sampled Parameters
%__________________________________________________________________________

% EDIT THE FOLLOWING VARIABLES, UNSAMPLED PARAMETERS, & ANY OTHER ARGS HERE

% Specify independent var(s) to pass to the model function
x = linspace(0,10,101);  % spatial domain for my testlinear function
% t = linspace(0,17,171);  % time domain for mypdes (commented example)

% Specify values of any unsampled parameters to pass to the function
unsampledps = 2; 

% Specify any additional solver parameters if needed: 
% w = 2; %spherical symmetry = 2 for pdepe.m to solve mypdes 

Simdata = struct; % initialize struct to store outputs for computing PRCCs

tic % start measuring time to solve equations to monitor progress

% Loop over the parameter sample pairs for Monte Carlo simulations
for j=1:N   
    sampledparams=A(j,1:M);
    fprintf('Parameters passed for sample: ');fprintf('%u',j);fprintf(' of ');fprintf('%u\n',N);
    
%  EDIT THE FOLLOWING FUNCTION CALL FOR YOUR OWN MODEL:
    Simdata(j).y=testlinear(x,sampledparams,unsampledps);

% % % EXAMPLE FOR COUPLED PDE MODEL SOLVED USING PDEPE:
% %     sol = pdepe(w,@mypdes,@mypdesic,@mypdesbc,x,t); 
% %     Simdata(j).c = sol(:,:,1); % solution to the equation for variable c
% %     Simdata(j).p = sol(:,:,2); % solution to the equation for variable p
% %     Simdata(j).r = sol(:,:,3); % solution to the equation for variable r

    toc
end

save(outFileName2, 'Simdata') % add Simdata to the saved .mat file


OutputOfInterest = zeros(N,length(x));
% OutputOfInterest = zeros(N,length(x)length(t))

for si = 1:N
% EDIT THE FOLLOWING TO SPECIFY OUTPUT DATA TO COMPARE:
    OutputOfInterest(si,:) = Simdata(si).y;
% % OutputOfInterest(si,:,:) =Simdata(si).r./(Simdata(si).c+Simdata(si).r)
end

% EDIT THE FOLLOWING STRING TO NAME OUTPUT DATA:
labelstring = 'y'; % consistent id for plot labels, filenames
% % labelstring = 'ratio'


%---Visualize the range of simulation results with errorbar plots---%

outPlot = plotSimulationOutput(x,N,OutputOfInterest,labelstring); 

pause(5)

figurelabel1=([outFileStr,'-N',num2str(N),'-',labelstring,'-ErrorbarPlot.fig']);
figurelabel2=([outFileStr,'-N',num2str(N),'-',labelstring,'-ErrorbarPlot.png']);
saveas(outPlot,figurelabel1);
saveas(outPlot,figurelabel2);

%%_________________________________________________________________________
%
% THIRD STEP: Compute & Plot Partial Rank Correlation Coefficients (PRCC)
%__________________________________________________________________________
    
tx_varied = input('Do you want to look at results that might change over time or space (''Y'' or ''N'')?: ');

if strcmp(tx_varied,'Y')

    % EDIT INDICES FOR EVALUATING RESULTS BELOW:

    % To look at variation in TIME: 

    % x_index = 5; % specify spatial index if a PDE model
    % x_index = []; % for time varied ODE model
%     [prcc]=TimeVariedPRCC(M,N,K,Simdata,t,x_index);

    % To look at variation in SPACE:

    % t_index= 10; % specify time index if a PDE model
    t_index =[]; % for a spatially varied ODE model
    prcc = VariedPRCC(M,N,A,OutputOfInterest,x,t_index);

    prccPlot = plotVariedPRCC(M,N,x,labelstring,parameters,prcc);
    
elseif strcmp(tx_varied,'N')

    % EDIT INDICES FOR EVALUATING RESULTS:
    x_index = 1; % if an ODE, leave the n/a index empty;
    t_index =[];  % if PDE, specify indices for both
    prcc = UnvariedPRCC(M,N,A,OutputOfInterest,t_index,x_index);
    
    paramNames={parameters.name};
    paramNames2=[parameters.name];

    prccPlot = plotUnvariedPRCC(M,N,labelstring,parameters,prcc);

end

figurelabel1=([outFileStr,'-N',num2str(N),'-', labelstring,'-PRCC.fig']); %assign an appropriate filename
figurelabel2=([outFileStr,'-N',num2str(N),'-', labelstring,'-PRCC.png']);
pause(5) %Time to dock/maximize the figure before it saves, if you prefer.
saveas(prccPlot,figurelabel1); %save fig
saveas(prccPlot,figurelabel2); %save png
