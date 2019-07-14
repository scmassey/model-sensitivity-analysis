function [prcc]=UnvariedPRCC(M,N,A,OutputOfInterest,t_index,x_index)

    K=zeros(N,M+1); % number of samples drawn N x number of model parameters (M) + 1 output column
    K([1:N],[1:M])=A([1:N],[1:M]); %Taking the information from A and putting it into K

	prcc=zeros(M+1,1); % will hold PRCC statistic for each parameter for chosen point (x=0, or otherwise)
    % studentT=zeros(M+1,1); % will hold the studentT statistic
    
    for n=1:N %loop through each of the N samples
        K(n,M+1)=OutputOfInterest(n,x_index); %Output results go in the added last column, the M+1st column.
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

    % for w=1:M %iterate thru the M parameters to calculate PRCCs for each with ratio
    %     prcc(w)=(-B(w,M+1))/sqrt(B(w,w)*B(M+1,M+1)); % the PRCC between the wth parameter and the ratio
    %     studentT(w)=prcc(w)*sqrt((N-2)/(1-prcc(w))); % the studentT statistic corresponding to the PRCC
    % end
    % %Each studentT is the distribution (showing the significance of) the corresponding PRCC value

end