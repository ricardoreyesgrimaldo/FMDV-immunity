function [opti,maxlikelihood,gradint,H,ci1,ci2,ci3,sigma]=mleestimator(initialcondition)
tic
% initialcondition is a 3-dimensional vector that provides 
datafilename='SAT.csv'; %Provide the name of data file
threshold=1.7;

data = xlsread(datafilename); %Obtain the data from the csv file.
b=find(~any(isnan(data(:,3)),2)); %Detect where there is missing data
data1=data(b,:); %Ignores where there is data missing
data1(:,3:5)=(data1(:,3:5)>=threshold)+1; %Verify data that is under 
                                          %threshold value
%sat=3;
for sat=1:3
    f=@(v) minusloglikelihood2(v(1),v(2),v(3),data1,sat);
    %Constraints for minimizer
    A=[-1,0,0;0,-1,0];
    b1=[0;0];
    Aeq=[];
    beq=[];
    lb=[0,0,0];
    ub=[Inf,Inf,Inf];
    %options=optimoptions(@fmincon,'Display','off');
    [opt,fval,exitflag,output,lambda,gradint,H]=fmincon(f,initialcondition(sat,:),A,b1,Aeq,beq,lb,ub);
    opti(sat,:)=opt;
    maxlikelihood(sat,:)=fval;
    sigma=inv(H)*(1/2);
    z=norminv(0.975);
    ci1(sat,:)=[opt(1)-z*sigma(1,1),opt(1)+z*sigma(1,1)];
    ci2(sat,:)=[opt(2)-z*sigma(2,2),opt(2)+z*sigma(2,2)];
    ci3(sat,:)=[opt(3)-z*sigma(3,3),opt(3)+z*sigma(3,3)];
end
toc
sprintf('mle estimator model 2 done')
end