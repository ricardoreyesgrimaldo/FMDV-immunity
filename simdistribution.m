%optis=[0.669900429045425,0.146556249607274,0.089120132317655;
%    0.273349126394630,0.383014310604545,0.206566501840294;
%    0.381970241054393,0.178230681553107,0.362475586639337];
%co=(0);
%[high1,count1,low1,count_low1]=distributionsat4(co);
function simdistribution(optimodel1,optimodel2,highmodel1,countmodel1,lowmodel1,count_lowmodel1,highmodel2,countmodel2,lowmodel2,count_lowmodel2,high1,count1,low1,count_low1)

for sat=1:3
    t=[1,3:18];
    alpha=optimodel2(sat,1);
    beta=optimodel2(sat,2);
    c=optimodel2(sat,3);
    a=optimodel1(sat,1);
    b=optimodel2(sat,2);
    for i=1:length(t)
        Ph1=ProbFMDV(t(i),a,b);
        p1(sat,i)=Ph1(2,2);
        %mu1=count1(i,sat)*p1(sat,i);
        %alpha11=mu1+1;
        %beta11=count1(i,sat)-mu1+1;
        lp1(sat,i)=binoinv(0.025,count1(i,sat),p1(sat,i))/count1(i,sat);
        %lp1(sat,i)=betainv(0.025,alpha11,beta11);
%         up1(sat,i)=betainv(0.975,alpha11,beta11);
        up1(sat,i)=binoinv(0.975,count1(i,sat),p1(sat,i))/count1(i,sat);
        q1(sat,i)=Ph1(1,2);
        mul1=count_low1(i,sat)*q1(sat,i);
        alphal1=mul1+1;
        betal1=count_low1(i,sat)-mul1+1;
        lq1(sat,i)=betainv(0.025,alphal1,betal1);
        uq1(sat,i)=betainv(0.975,alphal1,betal1);
        %%%%%%%%%%%%
        Ph2=Phigh(t(i),1,alpha,beta,c,1,18);
        p2(sat,i)=Ph2(2,2);
        mu2=count1(i,sat)*p2(sat,i);
        alpha12=mu2+1;
        beta12=count1(i,sat)-mu2+1;
        lp2(sat,i)=betainv(0.025,alpha12,beta12);
        up2(sat,i)=betainv(0.975,alpha12,beta12);
        q2(sat,i)=Ph2(1,2);
        mul2=count_low1(i,sat)*q2(sat,i);
        alphal2=mul2+1;
        betal2=count_low1(i,sat)-mul2+1;
        lq2(sat,i)=betainv(0.025,alphal2,betal2);
        uq2(sat,i)=betainv(0.975,alphal2,betal2);
    end
end
for sat=1:3
    figure()
    hold on
    plot(p1(sat,:))
    plot(p2(sat,:))
    plot(high1(:,sat)./count1(:,sat))
    curve1=lp1(sat,:);
    plot(lp1(sat,:))
    curve2=up1(sat,:);
    plot(up1(sat,:))
    x=[1:length(lp1(sat,:))];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'g', 'FaceAlpha',0.2);
    curve3=lp2(sat,:);
    plot(lp2(sat,:))
    curve4=up2(sat,:);
    plot(up2(sat,:))
    x=[1:length(lp2(sat,:))];
    x2 = [x, fliplr(x)];
    inBetween = [curve3, fliplr(curve4)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.2);
end
for sat=1:3
    figure()
    hold on
    plot(q1(sat,:))
    plot(q2(sat,:))
    plot(1-low1(:,sat)./count_low1(:,sat))
    curve1=lq1(sat,:);
    plot(lq1(sat,:))
    curve2=uq1(sat,:);
    plot(uq1(sat,:))
    x=[1:length(lq1(sat,:))];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, 'g', 'FaceAlpha',0.2);
    curve3=lq2(sat,:);
    plot(lq2(sat,:))
    curve4=uq2(sat,:);
    plot(uq2(sat,:))
    x=[1:length(lq2(sat,:))];
    x2 = [x, fliplr(x)];
    inBetween = [curve3, fliplr(curve4)];
    fill(x2, inBetween, 'r', 'FaceAlpha',0.2);
end
