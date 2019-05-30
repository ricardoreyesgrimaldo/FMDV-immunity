function [high,count,low,count_low]=distributionsat4(co)
%data = xlsread('data-jan.csv');
data = xlsread('Data\SAT.csv');
a=unique(data(:,1));%Number of animals
b=unique(data(:,2));%Capture Time
n=length(unique(data(:,2)));
count=co;

high=zeros(max(b),3);
low=zeros(max(b),3);
count=zeros(max(b),3);
count_low=zeros(max(b),3);
for k=1:length(a)
    A=data(:,1)==a(k);
    for m=1:3
    %B=data(A,2:5+m:7+m);
    B=data(A,2:m:2+m);
    B=B(~isnan(B(:,2)),:);
    if (size(B, 1) > 0)
        post=(B(:,2)>=1.7);
        for j=1:size(B,1)
            i=B(j,1)-B(1,1)+1;
            if post(1)==1 
                if post(j)==1
                    high(i,m)=high(i,m)+1;
                end
                count(i,m)=count(i,m)+1;
            else % post(1)==0
                if post(j)==0
                    low(i,m)=low(i,m)+1;
                end
                count_low(i,m)=count_low(i,m)+1;
            end
        end
    end
    end
end