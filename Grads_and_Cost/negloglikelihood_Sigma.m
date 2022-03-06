function res = negloglikelihood_Sigma(data,w,mu,Sigma,reparam)
    %Sigma:Kx1cell with DxDarray, mu:KxD
    N=size(data,1);
    K=size(Sigma,1);
    D=size(mu,2);
    if size(w,2)==(K-1)
        w=[w,0];
        p=exp(w-max(w)) ./ sum(exp(w-max(w)));
    else 
        p=w;
    end
    mvns = zeros(N,K);
    for i=1:K
        if any(isnan(Sigma{i}),"all") || any(isinf(Sigma{i}),"all")
            res=NaN;
            return
        else
            eigs=eig(Sigma{i});
            if any(eigs<=0)
                res=NaN;
                return
            end 
        end 
        mvns(:,i)=reshape(mvnpdf(data,mu(i,:),Sigma{i}),N,1);
    end
    if reparam
        factor=2*pi*exp(0.5);
    else
        factor=1;
    end 
    res=0;
    for i=1:N
        a=10e-200;
        for j=1:K
            a = a+p(j)*factor*mvns(i,j);
        end
        a=log(a);
        res=res-a;
    end
    