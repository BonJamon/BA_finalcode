function grad_riem=grad_a_2(a,mvns)
    K=size(a,1);
    grad_euclid=reshape(mvns,K,1);
    Z=grad_euclid.*a;
    alpha=ones(1,K)*Z;
    grad_riem=Z-(ones(K,1)*alpha).*a;
