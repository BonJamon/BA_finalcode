function [w,mu,X]=initialize_Parameters(data,K,seed)
    rng(seed)
    w=[randn(1,K-1),0];
    D=size(data,2);
    [idx,mu]=kmeans(data,K);
    Sigma=cell(K,1);
    for i=1:K
        Sigma{i}=cov(reshape(data(idx==i,:),[],D));
        [T,flag]=cholcov(Sigma{i},0);
        if flag>0
            Sigma{i}=eye(D);
        end 
    end
    X=get_A_from_S(Sigma);
    X={X};