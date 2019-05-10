function l=minusloglikelihood(a,b,data,sat) % Minus Log Likelihood Uses estimates of parameters a and b along with data(added sat as a parameter)

na=unique(data(:,1));%Number of animals
l=0;
for i=1:length(na)
%     A=data(:,1)==na(i);
%     B=data(A,2:3);
    A=data(:,1)==na(i);
    B=data(A,2:sat:2+sat);
    for j=2:size(B,1)
        t=B(j,1)-B(j-1,1);
        P=ProbFMDV(t,a,b);
        l=l+log(P(B(j-1,2),B(j,2)));
    end
end
l=-l;