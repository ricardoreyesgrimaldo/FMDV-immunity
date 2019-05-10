function P=ProbFMDV(t,a,b)
P=[1-(a/(a+b))*(1-exp(-t.*(a+b))),(a/(a+b))*(1-exp(-t.*(a+b)));
    (b/(a+b))*(1-exp(-t.*(a+b))),1-(b/(a+b))*(1-exp(-t.*(a+b)))];
end
