function res = EM(data,w,mu,X,N,K,D,maxiter)

    function gamma=E_step(data,p,mu,Sigma)
        gamma = zeros(N,K);
        mvns = zeros(N,K);
        for i=1:K
            mvns(:,i)=reshape(mvnpdf(data,mu(i,:),Sigma{i}),N,1);
        end
        for i=1:N
            sum_i=sum(p .* mvns(i,:));
            for j=1:K
                gamma(i,j)=p(j)*mvns(i,j)/sum_i;
            end
        end
    end

    function [p,mu,Sigma]=M_step(data,gamma)
        mu=zeros(K,D);
        Sigma=cellmat(K,1,D,D,0);
        
        n_k=zeros(1,K);
        for i=1:K
            n_k(i)=sum(gamma(:,i));
        end

        p = n_k / N;
        for j=1:K
            for i=1:N
                mu(j,:)=mu(j,:)+gamma(i,j)*data(i,:);
            end 
            mu(j,:)=mu(j,:)/n_k(j);
        end 

        
        for j=1:K
            for i=1:N
                vec=reshape(data(i,:)-mu(j,:),D,1);
                Sigma{j}=Sigma{j}+gamma(i,j)*(vec*vec');
            end 
            Sigma{j}=Sigma{j} / n_k(j);
            [T,flag]=cholcov(Sigma{j},0);
            if flag>0
                singular=true;
            end 
        end 
    end 
    
    %Start
    singular=false;
    p = exp(w-max(w)) ./sum(exp(w-max(w)));
    Sigma = get_S_from_A(X{1});
    costs = zeros(maxiter+1,1);
    times = zeros(maxiter+1,1);
    xs =struct("p",zeros(maxiter+1,K), "mu",zeros(maxiter+1,K,D),"Sigma",cell(maxiter+1,K,1));
    costs(1)=negloglikelihood_Sigma(data,p,mu,Sigma,false);
    xs(1)=struct("p",p,"mu",mu,"Sigma",{Sigma});
    for iter=2:maxiter+1
        timetic = tic();
        gamma=E_step(data,p,mu,Sigma);
        [p,mu,Sigma]=M_step(data,gamma);
        
        if singular
            xs=xs(1:iter-1);
            costs=costs(1:iter-1);
            times=times(1:iter-1);
            break
        end 
        
        times(iter)=times(iter-1)+toc(timetic);
        costs(iter)=negloglikelihood_Sigma(data,p,mu,Sigma,false);
        xs(iter)=struct("p",p,"mu",mu,"Sigma",{Sigma});
        if costs(iter-1)-costs(iter)<1e-6
            xs=xs(1:iter);
            costs=costs(1:iter);
            times=times(1:iter);
            break
        end
    end 
    res=struct("xs",xs,"costs",costs,"times",times,"singular",singular);         
end        