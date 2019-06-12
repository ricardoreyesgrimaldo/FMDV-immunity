clear; close all; home; format long g;
tic
load('mlemodel1.mat')
threshold=1.7;
% initialcondition is a 3-dimensional vector that provides 
datafilename='SAT.csv'; %Provide the name of data file
%threshold=1.7;

data = xlsread(datafilename); %Obtain the data from the csv file.
b=find(~any(isnan(data(:,3)),2)); %Detect where there is missing data
data1=data(b,:); %Ignores where there is data missing
data1(:,3:5)=(data1(:,3:5)>=threshold)+1; %Verify data that is under 
                                          %threshold value
%sat=1;
for sat=1:3
%% Initial data
% Parameters
nn      = 100;       % Number of samples for examine the AC
N       = 30000;     % Number of samples (iterations)
burnin  = 3000;      % Number of runs until the chain approaches stationarity
lag     = 10;        % Thinning or lag period: storing only every lag-th point
% Storage
theta   = zeros(2,N);      % Samples drawn from the Markov chain (States)
acc     = 0;               % Accepted samples

a1=opti(sat,1);
cia=ci1(sat,:);
b1=opti(sat,2);
cib=ci2(sat,:);

%%Distribution
Sigma=[.5;.05];
proposal_PDF = @(v,ci) prod(unifpdf(v,ci-0.5*Sigma,ci+0.5*Sigma));
sample_from_proposal_PDF = @(ci) unifrnd(ci-0.5*Sigma,ci+0.5*Sigma);
%p=@(v) -minusloglikelihood(v(1),v(2),data1,sat);
p=@(v) pr(v,data1,sat);
aa=eps;  bb=1;
tt=[a1;b1];
%% Marginals
X = (cia(1):(cia(2)-cia(1))/241:cia(2))';   nx = length(X);
Y = (cib(1):(cib(2)-cib(1))/241:cib(2))';   ny = length(Y);
[XX,YY] = meshgrid(X,Y);
pXY = @(x,y) exp(-minusloglikelihood(x',y',data1,sat));
pX = zeros(nx,1);
for i = 1:nx
   pX(i) = integral(@(y) pXY(repmat(X(i),1,length(y)),y), 0, 1,'ArrayValued',true);  % Marginal X
end
pY = zeros(ny,1);
for i = 1:ny
   pY(i) = integral(@(x) pXY(x,repmat(Y(i),1,length(x))), 0, 1,'ArrayValued',true);  % Marginal Y
end
%% M-H routine
for i = 1:burnin   % First make the burn-in stage
   [tt a] = MH_routine(tt,p,proposal_PDF,sample_from_proposal_PDF); 
end
for i = 1:N   % Cycle to the number of samples
   for j = 1:lag   % Cycle to make the thinning
      [tt a] = MH_routine(tt,p,proposal_PDF,sample_from_proposal_PDF);
   end
   theta(:,i) = tt;        % Store the chosen states
   acc        = acc + a;   % Accepted ?
end
accrate = acc/N;           % Acceptance rate
%% Autocorrelation
%nn = 100;          % Number of samples for examine the AC
pp = theta(1,1:nn);   pp2 = theta(1,end-nn:end);   % First ans Last nn samples in X
qq = theta(2,1:nn);   qq2 = theta(2,end-nn:end);   % First and Last nn samples in Y
% AC in X
[r lags]   = xcorr(pp-mean(pp), 'coeff');
[r2 lags2] = xcorr(pp2-mean(pp2), 'coeff');
% AC in Y
[r3 lags3] = xcorr(qq-mean(qq), 'coeff');
[r4 lags4] = xcorr(qq2-mean(qq2), 'coeff');

%% Plots
% Autocorrelation
figure;
subplot(2,2,1);   stem(lags, r);
title('Autocorrelation in X', 'FontSize', 14);
ylabel('AC (first 100 samples)', 'FontSize', 12);
subplot(2,2,3);   stem(lags2, r2);
ylabel('AC (last 100 samples)', 'FontSize', 12);
subplot(2,2,2);   stem(lags3, r3);
title('Autocorrelation in Y', 'FontSize', 14);
ylabel('AC (first 100 samples)', 'FontSize', 12);
subplot(2,2,4);   stem(lags4, r4);
ylabel('AC (last 100 samples)', 'FontSize', 12);
% Target function and samples 
for i=1:length(XX)
    for j=1:length(YY)
        Z(i,j) = p([XX(i,j) YY(i,j)]);
    end
end
Z = reshape(Z,length(YY),length(XX));
figure;
subplot(2,1,1);   % Target "PDF"
surf(X,Y,Z); grid on; shading interp;
xlabel('X', 'FontSize', 12);  ylabel('Y', 'FontSize', 12);  
title('f_{XY}(x,y)', 'FontSize', 12);
subplot(2,1,2);   % Distribution of samples
plot(theta(1,:),theta(2,:),'k.','LineWidth',1); hold on; 
contour(X,Y,Z,22,'--','LineWidth',2); colormap jet; axis tight;
xlabel('X', 'FontSize', 12);   ylabel('Y', 'FontSize', 12);
text(3,-3,sprintf('Acceptace rate = %g', accrate),'FontSize',11);
% Joint, Marginals and histograms
figure;
% Marginal Y
subplot(4,4,[1 5 9]);
[n2 x2] = hist(theta(2,:), ceil(sqrt(N))); 
barh(x2, n2/(N*(x2(2)-x2(1))));   hold on;               % Normalized histogram
set(gca,'XDir','reverse','YAxisLocation','right', 'Box','off'); 
xlabel('pY(y)', 'FontSize', 15) 
plot(pY/trapz(Y,pY),Y,'r-','LineWidth',2); axis tight;   % Normalized marginal
% Marginal X
subplot(4,4,[14 15 16]);
[n1 x1] = hist(theta(1,:), ceil(sqrt(N))); 
bar(x1, n1/(N*(x1(2)-x1(1))));   axis tight;   hold on;  % Normalized histogram
set(gca,'YDir','reverse','XAxisLocation','top','Box','off'); 
ylabel('pX(x)', 'FontSize', 15)
plot(X,pX/trapz(X,pX),'r-','LineWidth',2);               % Normalized marginal
% Distribution of samples
subplot(4,4,[2 3 4 6 7 8 10 11 12]);  
plot(theta(1,:),theta(2,:),'b.','LineWidth',1); axis tight; hold on; 
contour(X,Y,Z,22,'--','LineWidth',2); colormap summer; 
title('Distribution of samples', 'FontSize', 15);
xlabel('X', 'FontSize', 15);   ylabel('Y', 'FontSize', 15);
% Useful Matlab command
figure;
scatterhist(theta(1,:),theta(2,:),[ceil(sqrt(N)) ceil(sqrt(N))]); hold on; 
contour(X,Y,Z,22,'--','LineWidth',2); colormap summer;

save(['bayesianrun4_' num2str(sat) '.mat'])
end
%save('bayesianrun3.mat')
%%End
toc
%Use percentiles and median to obtain credible interval and MLE from
%Bayesian approach