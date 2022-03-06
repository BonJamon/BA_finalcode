function [w,mu,X]=initialize_random(K,D,seed)
    rng(seed)
    w=[randn(1,K-1),0];
    mu=randn(K,D);
    Sigma=cell(K,1);
    for i=1:K
        Sigma{i}=eye(D);
    end 
    X=get_A_from_S(Sigma);
    X={X};
end 