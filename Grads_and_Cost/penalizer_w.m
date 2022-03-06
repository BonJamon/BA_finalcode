function out = penalizer_w(w)
    K=size(w,2);
    c=1;
    out = c*sum(w)-K*c*log(sum(exp(w)));