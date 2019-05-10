function [opti,maxlikelihood]=mleestimatormodel1(threshold)
tic
% initialcondition is a 3-dimensional vector that provides 
datafilename='SAT.csv'; %Provide the name of data file
%threshold=1.7;

data = xlsread(datafilename); %Obtain the data from the csv file.
b=find(~any(isnan(data(:,3)),2)); %Detect where there is missing data
data1=data(b,:); %Ignores where there is data missing
data1(:,3:5)=(data1(:,3:5)>=threshold)+1; %Verify data that is under 
                                          %threshold value
%sat=3;
A=[2.822257, 2.012271 , 1.659891;0.543653, 1.191523 , 2.098627];
parfor sat=1:3
    f=@(v) minusloglikelihood(v(1),v(2),data1,sat);
    %options=optimoptions(@fmincon,'Display','off');
    [opt,fval]=fminunc(f,A(:,sat));
    opti(sat,:)=opt;
    maxlikelihood(sat,:)=fval;
%     sigma=inv(H)*(1/2);
%     z=norminv(0.975);
%     ci1(sat,:)=[opt(1)-z*sigma(1,1),opt(1)+z*sigma(1,1)];
%     ci2(sat,:)=[opt(2)-z*sigma(2,2),opt(2)+z*sigma(2,2)];
%     ci3(sat,:)=[opt(3)-z*sigma(3,3),opt(3)+z*sigma(3,3)];
end
toc
sprintf('mle estimator model 1 done')
end