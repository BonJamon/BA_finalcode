function S =get_S_from_mu_Sigma(mu,Sigma)
    K=size(Sigma,1);
    D=size(Sigma{1},1);
    S=cellmat(K,1,D+1,D+1);
    for i=1:K
        S{i}(1:D,1:D)=Sigma{i}+mu(i,:)'*mu(i,:);
        S{i}(1:D,D+1)=mu(i,:)';
        S{i}(D+1,1:D)=mu(i,:);
        S{i}(D+1,D+1)=1;
    end
    