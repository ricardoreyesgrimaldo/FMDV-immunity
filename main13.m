close all
clear all
format long
tic
runsize=5;
count=zeros(3,1);
countin=zeros(3,1);
counton=zeros(3,1);
count1=zeros(3,1);
countin1=zeros(3,1);
counton1=zeros(3,1);
na=50;
SATS=zeros(na*18,3*runsize+2);
%gx=zeros(3*runsize,101);
%gy=zeros(3*runsize,101);
data=zeros(na*18*runsize,5);
opti=zeros(runsize,3,3);
initialcondition=[0.669900360600708,0.146556927474911,0.089120169578574;
    0.273350090432041,0.383014797954857,0.206566767497954;
    0.381970188979518   0.178230531387283   0.362475512442551];
%initialcondition=[0.278370740957182,0.254473970566066,-1.081339737935759;-1.297009166123985,-0.959681369285760,-1.577133037826026;-0.962413248735407,-1.724673400288835,-1.014797518110715];
%[optis]=mleestimator(initialcondition);
optis=[0.669900429045425,0.146556249607274,0.089120132317655;0.273349126394630,0.383014310604545,0.206566501840294;0.381970241054393,0.178230681553107,0.362475586639337];
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
        SATS(:,1:2)=SAT(:,1:2);
        SATS(:,3*run+i-1)=SAT(:,i+2);
        %f=@(v) minusloglikelihood2(v(1),v(2),v(3),SAT,i);
%         Y=linspace(a(1)/2,100*a(1));
%         Z=linspace(b(1)/2,100*b(1));
%         W=zeros(length(Y),length(Z));
%         for j=1:length(Y)
%             for k=1:length(Z)
%                 W(j,k)=f([Y(j),Z(k)]);
%             end
%         end
%         contour(Y,Z,W)
%         colorbar
%         hold on
%         plot(a(1),b(1),'*')
%         figure()
%         pcolor(Y,Z,W)
%         colorbar
%         hold on 
%         plot(a(1),b(1),'*')

% %Constraints for minimizer
%     A=[-1,0,0;0,-1,0];
%     b1=[0;0];
%     Aeq=[];
%     beq=[];
%     lb=[0,0,0];
%     ub=[Inf,Inf,Inf];
%     %options=optimoptions(@fmincon,'Display','off');
%     [opt,cn,e,s,g,k,H]=fmincon(f,initialcondition,A,b1,Aeq,beq,lb,ub);

        
%         options2=optimoptions(@fsolve,'Display','off');
%         accur=50;
%         theta=0:pi/accur:2*pi;
%         r=zeros(1,length(theta));
%         for m=1:length(theta)
%             g=@(r) f([(exp(r))*cos(theta(m))+opt(1),(exp(r))*sin(theta(m))+opt(2)])-(f([opt(1),opt(2)])+2);
%             r(m)=exp(fsolve(g,-1,options2));  
%             if imag(r(m))== 0
%                 countm=countm+1;
%                 r1(countm)=r(m);
%                 theta1(countm)=theta(m);
%             end
%         end
%         gx(i+3*(run-1),:)=r1.*cos(theta1)+opt(1);
%         gy(i+3*(run-1),:)=r1.*sin(theta1)+opt(2);
% 
%          opti(run,i,:)=opt;
%         [in,on] = inpolygon(opt(1),opt(2),r1.*cos(theta1)+opt(1),r1.*sin(theta1)+opt(2));
%         count(i)=count(i)+double(in)+double(on);
%         countin(i)=countin(i)+double(in);
%         counton(i)=counton(i)+double(on);
%         [in1,on1] = inpolygon(a(i),b(i),r1.*cos(theta1)+opt(1),r1.*sin(theta1)+opt(2));
%         count1(i)=count1(i)+double(in1)+double(on1);
%         countin1(i)=countin1(i)+double(in1);
%         counton1(i)=counton1(i)+double(on1);
    end
end

% for i=1:3
%     figure()
%     plot(a(i),b(i),'r*','LineWidth',5)
%     hold on
%     for k=1:runsize
%         fill(gx(3*(k-1)+i,:),gy(3*(k-1)+i,:),rand(1,3),'FaceAlpha',0.1)
%     end
%     axis([0 2*a(i) 0 2*b(i)])
% end
countna=0;
for run=1:runsize
    for i=1:18*na
        data(18*na*(run-1)+i,1)=SATS(i,1)+(countna)*na;
        data(18*na*(run-1)+i,2)=SATS(i,2);
        data(18*na*(run-1)+i,3:5)=SATS(i,3*run:3*run+2);
    end
    countna=countna+1;
end
[high,count,low,count_low]=distributionsat3(data);
co=0;
[high1,count1,low1,count_low1]=distributionsat4(co);
h=high;
c=count;
t=[1,3:19];
alpha=zeros(1,length(t));
beta=zeros(1,length(t));
for j=1:3
%j=1;
    figure()
    t=[1,3:19];
    hold on
    plot(high(:,j)./count(:,j),'--','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    plot(high1(:,j)./count1(:,j),'r','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    title({'Proportion of animals that started with high antibody levels','and had high antibody level at capture $k$'},'Interpreter','Latex')
    xlabel('Capture number')
    ylabel({'Percentage with respect of total number of animals','with high antibody level at capture $k$'},'Interpreter','Latex')
    SAT1=zeros(1,length(t));
    p=zeros(2,length(t));
    for i=1:length(t)
        alpha(i)=h(i,j)+1;
        beta(i)=c(i,j)-high(i,j)+1;
        %p(i)=betapdf(t(i)/max(t),alpha(i),beta(i));
        p(1,i)=betainv(0.975,alpha(i),beta(i));
        p(2,i)=betainv(0.025,alpha(i),beta(i));
        SAT1(i)=binopdf(t(i),c(i,j),p(i));
    end
    %figure()
    curve1=p(1,:);
    plot(p(1,:))
    curve2=p(2,:);
    plot(p(2,:))
    x=[1:18];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'g', 'FaceAlpha',0.2);
end
for j=1:3
%j=1;
    figure()
    t=[1,3:19];
    hold on
    plot(1-low(:,j)./count_low(:,j),'--','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    plot(1-low1(:,j)./count_low1(:,j),'r','LineWidth',2)% SAT1 Proportion of individuals that started high and stayed high
    title({'Proportion of animals that started with low antibody levels','and had high antibody level at capture $k$'},'Interpreter','Latex')
xlabel('Capture number')
ylabel({'Percentage with respect of total number of animals','with high antibody level at capture $k$'},'Interpreter','Latex')
    SAT1=zeros(1,length(t));
    p=zeros(2,length(t));
    for i=1:length(t)
        alpha(i)=low(i,j)+1;
        beta(i)=count_low(i,j)-low(i,j)+1;
        %p(i)=betapdf(t(i)/max(t),alpha(i),beta(i));
        p(1,i)=betainv(0.975,alpha(i),beta(i));
        p(2,i)=betainv(0.025,alpha(i),beta(i));
        SAT1(i)=binopdf(t(i),c(i,j),p(i));
    end
    %figure()
    curve1=1-p(1,:);
    plot(1-p(1,:))
    curve2=1-p(2,:);
    plot(1-p(2,:))
    x=[1:18];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'g', 'FaceAlpha',0.1);
end

save('run5.mat')
toc