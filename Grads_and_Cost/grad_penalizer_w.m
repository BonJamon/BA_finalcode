function grad = grad_penalizer_w(w)
    c=1;
    K=size(w,2);
    T=sum(exp(w));
    grad=c*(1-K/T*(exp(w(1:K-1))));
    
    