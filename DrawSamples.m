function [parameters, A]=DrawSamples(M,N) 

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
        
        x1 = sqrt(2*ssigma)*erfinv(-0.9999)+mu;

        for n=1:N

            if n~=N
                x2 = sqrt(2*ssigma)*erfinv(2*n*intervalarea-1)+mu;
            elseif n==N
                x2 = sqrt(2*ssigma)*erfinv(0.9999)+mu;
            end

            parameters(m).sample(n)=rand(1)*abs(x2-x1)+x1;

            B(n,m)=parameters(m).sample(n); %store the sample value in the B matrix

            x1=x2;

        end


        fprintf(['All samples drawn for parameter ', parameters(m).name,'.\n']);
        
    elseif strcmp(parameters(m).distribution,'triangle')==1

        mmin  = input('Enter the min: '); %lower bound of triangle (the parameter distribution)
        mmax  = input('Enter the max: '); %upper bound of the parameter distribution
        mmode = input('Enter the mode: ');
        
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

        mmin = input('Enter minimum value of the uniform parameter: ');
        mmax = input('Enter maximum value of the uniform parameter: ');

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

end
