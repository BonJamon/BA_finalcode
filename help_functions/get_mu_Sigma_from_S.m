function [mu,Sigma] = get_mu_Sigma_from_S(S)
    K=size(S,1);
    D=size(S{1},1)-1;
    Sigma=cellmat(K,1,D,D);
    mu=zeros(K,D);
    for i=1:K
        m=S{i}(1:D,D+1)/S{i}(D+1,D+1);
        Sigma{i}=S{i}(1:D,1:D)-m*m'/S{i}(D+1,D+1);
        mu(i,:)=m;
    end
