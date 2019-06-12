function x=pr(v,data1,sat)
if any(v<0)==1
    x=-Inf;
else
    x=-minusloglikelihood(v(1),v(2),data1,sat);
end
end