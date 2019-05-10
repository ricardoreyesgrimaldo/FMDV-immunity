% Likelihood function 
function Ph=Phigh(t,t0,alpha,beta,c,tstart,tend)
a=(beta-alpha)/(tend-tstart);
b=((alpha-beta)/(tend-tstart))*tstart+alpha;
d=(a/2)*t0.^2+t0.*(b+c);
mu=exp(-(a/2)*t.^2-t.*(b+c)+d);%\mu^{-1} inverse integrating factor
if a==0
    Ph=ProbFMDV(t-t0,b,c);
else
    if a<0
%         gamma=c*sqrt(pi/(-2*a))*exp(-d-((b+c)^2/(2*a)));
%         nu=sqrt(-a/2)*(t+((b+c)/a));
%         nu0=sqrt(-a/2)*(t0+((b+c)/a));
%         p0=1+mu*(-1-gamma*(erf(nu)-erf(nu0)));
%         p1=1+mu*(1-1-gamma*(erf(nu)-erf(nu0)));
        gamma=c*sqrt(2*pi/(-a))*exp(-d-((b+c)^2/(2*a)));
        N=normcdf(t,-(b+c)/a,sqrt(-1/a));
        N0=normcdf(t0,-(b+c)/a,sqrt(-1/a));
        p0=1+mu*(-1-gamma*(N-N0));
        p1=1+mu*(-gamma*(N-N0));        
        Ph=[1-p0,p0;1-p1,p1];
    else
        gamma=c*sqrt(pi/(2*a))*exp(-d-((b+c)^2/(2*a)));
        nu=sqrt(a/2)*(t+((b+c)/a));
        nu0=sqrt(a/2)*(t0+((b+c)/a));
        N=erfi([nu;nu0]);
        %N0=erfi(nu0);
        N0=gamma*(N(1)-N(2));
        p0=1+mu*(-1-N0);%figure out what happens with erfi
        p1=1+mu*(-N0);
        Ph=[1-p0,p0;1-p1,p1];
    end
end
% Choose d to vanish exponent in gamma instead of cancelling mu??