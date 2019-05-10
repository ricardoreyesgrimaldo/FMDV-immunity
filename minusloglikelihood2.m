% Finds the minus log likelihood function for FMDV given specific data. 
% 

function l=minusloglikelihood2(alpha,beta,c,data,sat) 

na=unique(data(:,1));%Number of animals
time=unique(data(:,2));%Time captures
tstart=min(time);
tend=max(time);
l=0;
for i=1:length(na)
    A=data(:,1)==na(i);
    B=data(A,2:sat:2+sat);
    for j=2:size(B,1)
        t=B(j,1);
        t0=B(j-1,1);
        P=Phigh(t,t0,alpha,beta,c,tstart,tend);
        l=l+log(P(B(j-1,2),B(j,2)));
    end
end
l=-l;