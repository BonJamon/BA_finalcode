function out = lbfgs_manifold(data,w0,mu0,X0,N,K,D,penalized,maxiter,myrlbfgs)
%Extrem numerisch instabil für penalized version
    if penalized
        a=1;
        b=1;
        Delta=0.01*cov(data);
        k=0.01;
        lambd=mean(data);
        Psi=construct_Psi(a,b,Delta,k,lambd);
    end 
    singular=false;
    samecost=false;
    lsfailed=false;

    %Create Problem Structure
    symposM = sympositivedefinitefactory(D+1);
    symposM.transp = symposM.paralleltransp;
    %symposM.retr = symposM.exp;
    struct1 = struct('S', powermanifold(symposM,K), 'w', euclideanfactory(1,(K-1)));
    manifold = productmanifold(struct1);
    problem.M = manifold;

    function res = negloglikelihood(varstruct)  
        S = varstruct.S;%S sometimes not positive definite?
        w = varstruct.w;
        w = [w,0];

        N = size(data,1);
        y = [data, ones(N,1)];
        alpha= exp(w-max(w)) ./ sum(exp(w-max(w)));
        mu=zeros(K,D+1);
        res=negloglikelihood_Sigma(y,alpha,mu,S,true);
        if penalized
            for i=1:K
                res=res-penalizer_S(S{i},Psi);
            end 
            res=res-penalizer_w(w);
        end
    end

    function grad_struct = grad_negloglikelihood(varstruct)
        S = varstruct.S;
        w = varstruct.w;
        w = [w,0];
        y = [data, ones(N,1)];
  
        S_grad = cellmat(K,1,D+1,D+1);
        w_grad = zeros(1,K-1);
     
        p = exp(w-max(w)) ./ sum(exp(w-max(w)));

        mu=zeros(1,D+1);

        mvns = zeros(N,K);
        for i=1:K
            mvns(:,i)=reshape(mvnpdf(y,mu,S{i}),N,1);
        end
        for i=1:N
            sum_j=10e-200;
            for j=1:K
                sum_j = sum_j +p(j)*mvns(i,j);
            end
            for j=1:K
                S_grad{j} = S_grad{j} - (p(j)*grad_S(y(i,:),S{j},mvns(i,j))/sum_j);
            end 
            for j=1:(K-1)
                w_grad(j) = w_grad(j) - p(j)*(mvns(i,j)-sum(p.*mvns(i,:)))/ sum_j;
            end
                
        end 
        if penalized
            for i=1:K
                S_grad{i}=S_grad{i}-grad_penalizer_S(S{i},Psi);
            end
            w_grad = w_grad -grad_penalizer_w(w);
        end
        grad_struct = struct('S',{S_grad}, 'w', w_grad);     
    end
    
    % Define the problem cost function and its gradient.
    problem.cost  = @negloglikelihood;
    problem.grad = @grad_negloglikelihood;

    %Solve.
    function stats = true_negloglikelihood(problem,x,stats)
        [mu,Sigma]=get_mu_Sigma_from_S(x.S);
        stats.true_cost = negloglikelihood_Sigma(data,x.w,mu,Sigma,false);
        stats.x = x;
        stats.grad=grad_negloglikelihood(x);
        ls = stats.linesearch;
        if class(ls)=="double"
            stats.costevals = 0;
        else
            stats.costevals =ls.costevals+1;
        end 
    end
    %Reason for 3 attempts: When ls fails, minstepsize might still be
    %triggered, just sometimes after second fail
    attempts=0;
    function stopnow = stopcrit(problem,x,info,last)
        if last>1
            if info(last).cost==info(last-1).cost
                if attempts>2
                    samecost=true;
                else 
                    attempts=attempts+1;
                end
            else 
                attempts=0;
            end 
        end 
        stopnow = samecost||singular;
    end 
    w0 = w0(1:(K-1));
    Sigma0 = get_S_from_A(X0{1});
    S0 = get_S_from_mu_Sigma(mu0,Sigma0);
    %S0 = cellmat(K,1);
    %for n=1:K
    %    S0{n}=eye(D+1);
    %end 
    x0 = struct("S",{S0},"w",w0);
    options=struct("maxiter", maxiter,"statsfun",@true_negloglikelihood,...
        "minstepsize",1e-12,"memory",10,"verbosity",1,"tolgradnorm",1e-6,"stopfun",@stopcrit);
    if myrlbfgs
    	[x, xcost, info,options,singular,lsfailed] = my_rlbfgs(problem,x0,options); 
    else
        [x, xcost, info] = rlbfgs(problem,x0,options); 
    end 
    costs = [info.true_cost];
    times = [info.time];
    xs =[info.x];
    grads=[info.grad];
    stepsizes=[info.stepsize];
    costevals=[info.costevals];
    out=struct("xs", xs,"costs",costs,"times",times,"grads",grads,"stepsizes",...
        stepsizes,"gradnorms",[info.gradnorm],"costevals",costevals,"singular",...
        singular,"samecost",samecost,"lsfailed",lsfailed);
    
end