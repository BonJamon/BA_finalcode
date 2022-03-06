function out = lbfgs_reparam(data,w0,mu0,X0,N,K,D,penalized,maxiter,myrlbfgs)
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
    struct1 = struct('A', powermanifold(euclideanfactory(D+1,D+1),K), 'w', euclideanfactory(1,K-1));
    manifold = productmanifold(struct1);
    problem.M = manifold;

    function res = negloglikelihood(varstruct)  
        A = varstruct.A;
        w = varstruct.w;
        w = [w,0];
        S=get_S_from_A(A);

        N = size(data,1);
        y = [data, ones(N,1)];
        alpha= exp(w-max(w)) ./ sum(exp(w-max(w)));
        mu=zeros(K,D+1);
        res=negloglikelihood_Sigma(y,alpha,mu,S,true);
        %{
        logdetS = logdet(S);
        loglikevec = logsumexp(alpha,log_q(y,S,logdetS));

        res = -sum(loglikevec);
        %}
        if penalized
            for i=1:K
                res=res-penalizer_S(S{i},Psi);
            end 
            res=res-penalizer_w(w);
        end

    end

    function grad_struct = grad_negloglikelihood(varstruct)
        A = varstruct.A;
        w = varstruct.w;
        w = [w,0];
        y = [data, ones(N,1)];
  
        A_grad = cellmat(K,1,D+1,D+1);
        w_grad = zeros(1,K-1);
     
        p = exp(w-max(w)) ./ sum(exp(w-max(w)));

        mu=zeros(1,D+1);

        mvns = zeros(N,K);
        for i=1:K
            mvns(:,i)=reshape(mvnpdf(y,mu,A{i}*A{i}'),N,1);
        end
 
        for i=1:N
            sum_j=10e-200;
            for j=1:K
                sum_j = sum_j + p(j)*mvns(i,j);
            end
            for j=1:K
                A_grad{j} = A_grad{j} - (p(j)*grad_A(y(i,:),A{j},mvns(i,j))/sum_j);
            end 
            for j=1:(K-1)
                w_grad(j) = w_grad(j) - p(j)*(mvns(i,j)-sum(p.*mvns(i,:))) / sum_j;
            end
        end 
        if penalized
            for i=1:K
                A_grad{i}=A_grad{i}-grad_penalizer_A(A{i},Psi);
            end
            w_grad = w_grad -grad_penalizer_w(w);
        end
        grad_struct = struct('A',{A_grad}, 'w', w_grad);     
    end
    
    % Define the problem cost function and its gradient.
    problem.cost  = @negloglikelihood;
    problem.grad = @grad_negloglikelihood;
    
    %Solve.
    function stats = true_negloglikelihood(problem,x,stats)
        S=get_S_from_A(x.A);
        [mu,Sigma]=get_mu_Sigma_from_S(S);
        stats.true_cost = negloglikelihood_Sigma(data,x.w,mu,Sigma,false);
        stats.x=x;
        stats.grad=grad_negloglikelihood(x);
        ls = stats.linesearch;
        if class(ls)=="double"
            stats.costevals = 0;
        else
            stats.costevals =ls.costevals+1;
        end 
    end
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
    A0 = get_A_from_S(S0);
    x0=struct("A",{A0},"w",w0);
    options=struct("maxiter",maxiter,"statsfun",@true_negloglikelihood,...
        "minstepsize",1e-12,"memory",10,"verbosity",1,"tolgradnorm",1e-6,"stopfun",@stopcrit);
    if myrlbfgs
    	[x, xcost, info,options,singular,lsfailed] = my_rlbfgs(problem,x0,options); 
    else
        [x, xcost, info] = rlbfgs(problem,x0,options); 
    end 
    costs = [info.true_cost];
    times = [info.time];
    grads = [info.grad];
    stepsizes = [info.stepsize];
    costevals = [info.costevals];
    out = struct("xs",[info.x],"costs", costs, "times",times,"grads",grads,...
        "stepsizes",stepsizes,"gradnorms",[info.gradnorm],"costevals",costevals,...
        "singular",singular,"samecost",samecost,"lsfailed",lsfailed);
    
end