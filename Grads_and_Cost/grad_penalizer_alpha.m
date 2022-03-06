function grad_riem = grad_penalizer_alpha(a)
    %c=1 --> a^(c-1)=1
    K=size(a,1);
    grad_euclid=zeros(K,1);
    for i=1:K
        vec = a;
        vec(i)=1;
        grad_euclid(i)=prod(vec);
    end
    Z=grad_euclid.*a;
    alpha=ones(1,K)*Z;
    grad_riem=Z-(ones(K,1)*alpha).*a;