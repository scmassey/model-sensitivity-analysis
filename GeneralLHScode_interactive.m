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

% Initializing matrices to store samples, struct to store parameter names, distribution type, and sample values:
A=zeros(N,M); 
B=zeros(N,M);
parameters=struct;

intervalarea=1/N; % area of an interval = change in y for the interval when looking at the CDF

for m=1:M
    if mod(m,10)==1 && mod(m,100)~=11
        parameters(m).name = input(['Please specify the name of the ',num2str(m),'st parameter: ']);
    elseif mod(m,10)==2 && mod(m,100)~=12
        parameters(m).name = input(['Please specify the name of the ',num2str(m),'nd parameter: ']);
    elseif m==M
        parameters(m).name = input(['Please specify the name of the ',num2str(m),'th (last) parameter: ']);
    else
        parameters(m).name = input(['Please specify the name of the ',num2str(m),'th parameter: ']);
    end
    
    while ~ischar(parameters(m).name)
        fprintf('Not a character array. ');
        parameters(m).name = input('Please re-enter the parameter name as type char: ');
    end

    parameters(m).distribution = input(['Type of distribution for parameter ',parameters(m).name,' (''normal'',''triangle'', or ''uniform'' are supported): ']);  %If they all have different distributions, maybe draw the samples for each parameter at THIS step, vs the loop below..
    while ~strcmp(parameters(m).distribution,'normal')&&~strcmp(parameters(m).distribution,'triangle')&&~strcmp(parameters(m).distribution,'uniform')
        fprintf('Not a supported distribution type. ')
        parameters(m).distribution = input('Please enter ''normal'',''triangle'', or ''uniform'': ');
    end

    if strcmp(parameters(m).distribution,'normal')==1

        mu = input('Enter the mean: '); % TO DO: add check that these are type double - and similarly for other distbn inputs below
        ssigma = input('Enter the variance (square of standard deviation): '); 
        
        x1 = ssigma*sqrt(2)*erfinv(-0.9999)+mu;

        for n=1:N

            x2=ssigma*sqrt(2)*erfinv(2*intervalarea+erf(x1-mu)/(ssigma*sqrt(2)))+mu;

            parameters(m).sample(n)=rand(1)*abs(x2-x1)+x1;

            B(n,m)=parameters(m).sample(n); %store the sample value in the B matrix

            x1=x2;

        end

        fprintf(['All samples drawn for parameter ', parameters(m).name,'.\n']);
        
    elseif strcmp(parameters(m).distribution,'triangle')==1

        mmode = input('Enter the mode: ');
        mmin  = input('Enter the min: '); %lower bound of triangle (the parameter distribution)
        mmax  = input('Enter the max: '); %upper bound of the parameter distribution
        
        y1=0;
        for n=1:N      
            % Recall that the PDF and CDF for a triangle distribution are piecewise
            % defined, with parameters max, min, and mode (the peak of the triangle).
            
            % The area of the interval tells us the change in y we
            % want to step through the CDF to determine each intervals
            % change in x (width) since these will be narrower or wider
            % depending on proximity to the mode. 
            

            y1=intervalarea*(n-1); % lower bound of area segment being considered } looking at cdf read off area as y's, and convert to get x's
            y2=intervalarea*(n);   % upper bound of area segment being considered }       x's in cdf are the same as the x's in pdf! :)
        
            % Find the lower x-value of the sampling interval, x1:
            if y1<=((mmode-mmin)/(mmax-mmin)) % if y1 <= (mode-min)/(max-min), then min <= x1 <= mode
                x1=sqrt(y1*(mmax-mmin)*(mmode-mmin))+mmin;    

            elseif (1-((mmax-mmode)/(mmax-mmin)))<y1 % if y1 > (mode-min)/(max-min), then mode <= x1 <= max
                x1=mmax-sqrt((1-y1)*(mmax-mmin)*(mmax-mmode));   

            end

            % Find the upper x-value of the sampling interval:
            if y2<=((mmode-mmin)/(mmax-mmin)) % if y2 <= (mode-min)/(max-min), then min <= x2 <= mode
                x2=sqrt(y2*(mmax-mmin)*(mmode-mmin))+mmin; 

            elseif (1-((mmax-mmode)/(mmax-mmin)))<y2 % if y2 > (mode-min)/(max-min), then mode <= x2 <= max
                x2=mmax-sqrt((1-y2)*(mmax-mmin)*(mmax-mmode));  

            end

            % Draw a random sample from the interval:
            parameters(m).sample(n)=rand(1)*abs(x1-x2)+x1; %obtain the abs(difference between interval limits)and scale by rand(1) to find the "amount" to add to x1 (the lower limit of the interval)
                                            % i.e., if abs(x1-x2)=5, x1=33, x2=37,then we take rand(1) 
                                            % and get some fraction, say 1/4=.25, so that 5*.25=1.25,
                                            % and the sample = x1+1.25=34.25, falling within the interval,
                                            % as desired. 
                                            
            B(n,m)=parameters(m).sample(n); %store the sample value in the B matrix
        end
        fprintf(['All samples drawn for parameter ', parameters(m).name,'.\n']);
        
    elseif strcmp(parameters(m).distribution,'uniform')==1

        mmax = input('Enter maximum value of the uniform parameter: ');
        mmin = input('Enter minimum value of the uniform parameter: ');

        intervalwidth=(mmax-mmin)/N; 

        for n=1:N

            x1=mmin+intervalwidth*(n-1); % lower bound of sampling interval
            x2=mmin+intervalwidth*(n);   % upper bound of sampling interval

            parameters(m).sample(n)=rand(1)*abs(x2-x1)+x1;

            B(n,m)=parameters(m).sample(n); %store the sample value in the B matrix

        end

        fprintf(['All samples drawn for parameter ', parameters(m).name,'.\n']);

    end
end


% The values stored in the rows of B represent samples from the diagonal 
% of the hypercube; permute columns in order to get random pairings:
for m=1:M
    r=randperm(N);     %randomly rearrange the columns
    for n=1:N
        s=r(n);
        A(n,m)=B(s,m); %store new arrangement; rows of A are randomly paired samples
    end
end

%__________________________________________________________________________
%
% SECOND STEP: Solve Model Function with Sampled Parameters
%__________________________________________________________________________

% TO DO: Need user inputs here to note the model function and args and add checks
% myfun =  input('Please specify model function: '); % THIS ISN'T WORKING - NOT SURE HOW TO DO THIS...

% EDIT THE FOLLOWING VARIABLES, UNSAMPLED PARAMETERS, AND ANY OTHER ARGS  
% YOU PLAN TO PASS TO YOUR FUNCTION HERE IN THIS FILE--*IF* YOU PREFER THAT 
% TO COMMAND LINE INPUT: 
x = input('Please specify the independent variables of the function (e.g. x): ');
unsampledps = input('Please specify values of any unsampled parameters to pass to the function in a vector (can be empty): ');

% TO DO: Replace with input args? (see above lines)
% %Initialize space and time for pdepe:
% %these are some that I've used for my model; you can enter whatever you
% %need for the solver you are using.
% w = 2; %spherical symmetry = 2
% x = linspace(0,1,200); % 1 cm spherically symmetric domain 
% t0=0; % start time 0
% tf=17; % end time 17 days
% t = linspace(t0,tf,171);% time vector w/ number of time steps desired
% ninitx=.03; %cm radius of initial sphere of influence of injection
% nt=t;
% nx=x;

Simdata = struct; %place where all output data will be stored for computing PRCCs 

tic % start measuring time to solve equations
for j=1:N   
    sampledparams=A(j,1:M);
    fprintf('Parameters passed for sample: ');fprintf('%u',j);fprintf(' of ');fprintf('%u\n',N);
    
%  EDIT THE FOLLOWING FUNCTION CALL (can name multiple outputs as in commented example below
%  Also make sure everything passes to your function correctly and all initial conditions,
%  etc. are specified before entering this for loop over j.
    Simdata(j).y=testlinear(x,sampledparams,unsampledps);
% %     sol = pdepe(w,@kineticpdgfpdes,@pdgfpdesic,@pdgfpdesbc,nx,nt); %make sure the right set of pdes are used
% %     Simdata(j).c = sol(:,:,1); % assign these appropriately to store the data from your solution vector in the struct
% %     Simdata(j).p = sol(:,:,2); % For example, the part of the solution corresponding to the c pde I called Simdata(j).c
% %     Simdata(j).r = sol(:,:,3); % and the (j) index is required to hold each sample. To call a specific point:
% %     Simdata(j).n = sol(:,:,4); % Simdata(j=3).c(timestep=25,spacestep=50), or whatever makes sense for your equations.
    toc
end
% TO DO: have user specify file name to save - but also want to do this
% earlier, before start running code, and maybe even save updates every 
% several iterations?

% EDIT THE FOLLOWING FILE NAME FOR SAVING YOUR SIMULATION DATA.
save LHS_testlinear %give the workspace an appropriate unique name

%__________________________________________________________________________
%
% THIRD STEP: Compute the Partial Rank Correlation Coefficients (PRCC)
%__________________________________________________________________________

K=zeros(N,M+1); % number of samples drawn N x number of model parameters (M) + 1 output column
R=zeros(N,M+1); % number of samples drawn N x number of model parameters (M) + 1 output column

K([1:N],[1:M])=A([1:N],[1:M]); %Taking the information from A and putting it into K

labelstring = input('What would you like to call model output in plots and filenames? (please enter as char): ');
while ~ischar(labelstring)
    fprintf('Not a character array. ');
    labelstring = input('Please re-enter the output label as type char: ');
end
    
    % results_varied = input('Do you want to look at results that might change over time or space (''Y'' or ''N'')?: ');

    % if strcmp(results_varied,'Y') && length(Simdata.y(:,1))>1 && length(Simdata.y(1,:))>1
    %     t_or_x= input('Note that analyzing both are unsupported, please choose...')


time_varied = input('Do you want to look at results that might change over time (''Y'' or ''N'')?: ');

% space_varied = input('Are your results spatially varied (''Y'' or ''N'')?: ');
% 
% if strmp(space_varied,'Y')
%     space = input('Please enter location for analysis :');
% end

if strcmp(time_varied,'Y')
    prcc=zeros(M+1,length(days)); % will hold PRCC statistic for each parameter for each time at chosen spatial point (x=0, or otherwise)
    studentT=zeros(M+1,length(days)); % will hold the studentT statistic
    
    space_varied = input('Are your results also spatially varied (''Y'' or ''N'')?: ');
    
    if strcmp(space_varied,'Y')
        space = input('Please enter location for analysis :');
    end
    
    % EDIT THE FOLLOWING - steps loop works for one independent variable, 
    % not two. Adjust as needed if you want to specify a particular
    % dimension/independent var.
    for steps=1:length(x)  % loop through time (or space) - TO DO: would like to be able to do both and get surfaces..

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


        for w=1:M %iterate thru the M parameters to calculate PRCCs for each with ratio
            prcc(w,steps)=(-B(w,M+1))/sqrt(B(w,w)*B(M+1,M+1)); % the PRCC between the wth parameter and the ratio
            studentT(w,steps)=prcc(w,steps)*sqrt((N-2)/(1-prcc(w,steps))); % the studentT statistic corresponding to the PRCC
        end
        %Each studentT is the distribution (showing the significance of) the
        %corresponding gamma/PRCC. 

    end %end the loop over steps
    
elseif strcmp(time_varied,'N')
    spot = input('Specify index at which to analyze output: ');
    while rem(M,1)~=0||M<=0
        spot = input('Index should be a positive integer. Please re-enter the index: ');
    end
    
    prcc=zeros(M+1,1); % will hold PRCC statistic for each parameter for chosen point (x=0, or otherwise)
    studentT=zeros(M+1,1); % will hold the studentT statistic
    
    for n=1:N %loop through each of the N samples
        output=Simdata(n).y(spot);
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