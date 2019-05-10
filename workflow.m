%Likelihood 1 uses parameters (a,b)
%Likelihood 2 uses parameters (alpha,beta,c)
close all
clear all
format long
tic
runsize=10;
na=50;
threshold=1.7;
SATSmodel1=zeros(na*18,3*runsize+2); %Data for Likelihood 1
SATSmodel2=zeros(na*18,3*runsize+2); %Data for Likelihood 2
[optimodel1,maxlikelihoodmodel1]=mleestimatormodel1(threshold);
initialcondition=[0.669900408805846,0.146555999855608,0.089120104390880;0.273350391523738,0.383014580896015,0.206566766121656;0.381970288754607,0.178230376898277,0.362475461538944];
%initialcondition=[0.669900210684230,0.146556914599346,0.089120152600158;0.273348579589497,0.383014244202118,0.206566543917455;0.381970585316840,0.178230490886034,0.362475611573978];
%initialcondition=[0.669900429045425,0.146556249607274,0.089120132317655;0.273349126394630,0.383014310604545,0.206566501840294;0.381970241054393,0.178230681553107,0.362475586639337];
[optimodel2,maxlikelihoodmodel2]=mleestimator(initialcondition);
a=(optimodel1(:,1))';
b=(optimodel1(:,2))';
for run=1:runsize
    SAT=zeros(na*18,5);
%     a=[2.822257, 2.012271 , 1.659891];
%     b=[0.543653, 1.191523 , 2.098627];
    for i=1:3 %Generates 50 animals randomly, by using estimation from MLE for all the SAT's 
    %i=1;
        countm=0;
        X=zeros(na,18);
        t=1:18;
        P=ProbFMDV(1,a(i),b(i));
        for N=1:na
            X(N,1)=binornd(1,1/2);
            for k=1:17
                p=P(X(N,k)+1,2);
                X(N,k+1)=binornd(1,p);
            end
        end
        c=1;
        for N=1:size(X,1)
            for k=1:size(X,2)
                SAT(c,1)=N;
                SAT(c,2)=k;
                SAT(c,i+2)=X(N,k)+1;
                c=c+1; 
            end
        end
        SATSmodel1(:,1:2)=SAT(:,1:2);
        SATSmodel1(:,3*run+i-1)=SAT(:,i+2);
    end
end
optis=optimodel2;
for run=1:runsize
    SAT=zeros(na*18,5);
    for i=1:3 %Generates 50 animals randomly, by using estimation from MLE for all the SAT's 
    %i=1;
        countm=0;
        X=zeros(na,18);
        t=1:18;
        tstart=min(t);
        tend=max(t);
        for N=1:na
            X(N,1)=binornd(1,1/2);
            for k=1:17
                P=Phigh(t(k+1),t(k),optis(i,1),optis(i,2),optis(i,3),tstart,tend);%%%%%%
                p=P(X(N,k)+1,2);
                X(N,k+1)=binornd(1,p);
            end
        end
        c=1;
        for N=1:size(X,1)
            for k=1:size(X,2)
                SAT(c,1)=N;
                SAT(c,2)=k;
                SAT(c,i+2)=X(N,k)+1;
                c=c+1; 
            end
        end
        SATSmodel2(:,1:2)=SAT(:,1:2);
        SATSmodel2(:,3*run+i-1)=SAT(:,i+2);
        
    end
end
countna=0;
for run=1:runsize
    for i=1:18*na
        datamodel1(18*na*(run-1)+i,1)=SATSmodel1(i,1)+(countna)*na;
        datamodel1(18*na*(run-1)+i,2)=SATSmodel1(i,2);
        datamodel1(18*na*(run-1)+i,3:5)=SATSmodel1(i,3*run:3*run+2);
        datamodel2(18*na*(run-1)+i,1)=SATSmodel2(i,1)+(countna)*na;
        datamodel2(18*na*(run-1)+i,2)=SATSmodel2(i,2);
        datamodel2(18*na*(run-1)+i,3:5)=SATSmodel2(i,3*run:3*run+2);
    end
    countna=countna+1;
end
[highmodel1,countmodel1,lowmodel1,count_lowmodel1]=distributionsat3(datamodel1);
[highmodel2,countmodel2,lowmodel2,count_lowmodel2]=distributionsat3(datamodel2);
co=0;
[high1,count1,low1,count_low1]=distributionsat4(co);
h1=highmodel1;
c1=countmodel1;
h2=highmodel2;
c2=countmodel2;
t=[1,3:19];
alpha1=zeros(1,length(t));
beta1=zeros(1,length(t));
alpha2=zeros(1,length(t));
beta2=zeros(1,length(t));
for j=1:3
%j=1;
    figure()
    t=[1,3:19];
    hold on
    plot(highmodel1(:,j)./countmodel1(:,j),'--','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    plot(highmodel2(:,j)./countmodel2(:,j),'+','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    plot(high1(:,j)./count1(:,j),'r','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    title({'Proportion of animals that started with high antibody levels','and had high antibody level at capture $k$'},'Interpreter','Latex')
    xlabel('Capture number')
    ylabel({'Percentage with respect of total number of animals','with high antibody level at capture $k$'},'Interpreter','Latex')
    %SAT1=zeros(1,length(t));
    p1=zeros(2,length(t));
    p2=zeros(2,length(t));
    for i=1:length(t)
        alpha1(i)=h1(i,j)+1;
        beta1(i)=c1(i,j)-highmodel1(i,j)+1;
        alpha2(i)=h2(i,j)+1;
        beta2(i)=c2(i,j)-highmodel2(i,j)+1;
        %p(i)=betapdf(t(i)/max(t),alpha(i),beta(i));
        p1(1,i)=betainv(0.975,alpha1(i),beta1(i));
        p1(2,i)=betainv(0.025,alpha1(i),beta1(i));
        p2(1,i)=betainv(0.975,alpha2(i),beta2(i));
        p2(2,i)=betainv(0.025,alpha2(i),beta2(i));
        %SAT1(i)=binopdf(t(i),c(i,j),p(i));
    end
    %figure()
    curve1=p1(1,:);
    plot(p1(1,:))
    curve2=p1(2,:);
    plot(p1(2,:))
    x=[1:18];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'g', 'FaceAlpha',0.2);
    curve3=p2(1,:);
    plot(p2(1,:))
    curve4=p2(2,:);
    plot(p2(2,:))
    x=[1:18];
    x2 = [x, fliplr(x)];
    inBetween = [curve3, fliplr(curve4)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.2);
end
for j=1:3
%j=1;
    figure()
    t=[1,3:19];
    hold on
    plot(1-lowmodel1(:,j)./count_lowmodel1(:,j),'--','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    plot(1-lowmodel2(:,j)./count_lowmodel2(:,j),'+','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    plot(1-low1(:,j)./count_low1(:,j),'r','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    title({'Proportion of animals that started with low antibody levels','and had high antibody level at capture $k$'},'Interpreter','Latex')
xlabel('Capture number')
ylabel({'Percentage with respect of total number of animals','with high antibody level at capture $k$'},'Interpreter','Latex')
    p1=zeros(2,length(t));
    p2=zeros(2,length(t));
    for i=1:length(t)
        alpha1(i)=lowmodel1(i,j)+1;
        beta1(i)=count_lowmodel1(i,j)-lowmodel1(i,j)+1;
        alpha2(i)=lowmodel2(i,j)+1;
        beta2(i)=count_lowmodel2(i,j)-lowmodel2(i,j)+1;
        %p(i)=betapdf(t(i)/max(t),alpha(i),beta(i));
        p1(1,i)=betainv(0.975,alpha1(i),beta1(i));
        p1(2,i)=betainv(0.025,alpha1(i),beta1(i));
        p2(1,i)=betainv(0.975,alpha2(i),beta2(i));
        p2(2,i)=betainv(0.025,alpha2(i),beta2(i));
        
    end
    %figure()
    curve1=1-p1(1,:);
    plot(1-p1(1,:))
    curve2=1-p1(2,:);
    plot(1-p1(2,:))
    x=[1:18];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'g', 'FaceAlpha',0.2);
    curve3=1-p2(1,:);
    plot(1-p2(1,:))
    curve4=1-p2(2,:);
    plot(1-p2(2,:))
    x=[1:18];
    x2 = [x, fliplr(x)];
    inBetween = [curve3, fliplr(curve4)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.2);
end

Likelihoodtest=maxlikelihoodmodel1./maxlikelihoodmodel2;
simdistribution(optimodel1,optimodel2,highmodel1,countmodel1,lowmodel1,count_lowmodel1,highmodel2,countmodel2,lowmodel2,count_lowmodel2,high1,count1,low1,count_low1)
save('run6.mat')
Likelihoodtest
toc
sprintf('workflow done')