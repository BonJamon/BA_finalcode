function out = lbfgs_a(data,w0,mu0,X0,N,K,D,maxiter,myrlbfgs)
    %Create Problem Structure
    multM = multinomialfactory(K,1);
    multM.retr = multM.exp;
    struct1 = struct('X', powermanifold(euclideanfactory(D,D),K),'mu',euclideanfactory(K,D), 'a', multM);
    manifold = productmanifold(struct1);
    problem.M = manifold;
    singular=false;
    samecost=false;
    lsfailed=false;

    function res = negloglikelihood(varstruct)  
        X = varstruct.X;
        mu = varstruct.mu;
        p = reshape(varstruct.a,1,K);
        Sigma=get_S_from_A(X);
        res=negloglikelihood_Sigma(data,p,mu,Sigma,false);
    end

    function grad_struct = grad_negloglikelihood(varstruct)
        X = varstruct.X;
        Sigma = get_S_from_A(X);
        mu = varstruct.mu;
        p = varstruct.a;

        X_grad = cellmat(K,1,D,D);
        mu_grad = zeros(K,D);
        a_grad = zeros(K,1);

        mvns = zeros(N,K);
        for i=1:K
            mvns(:,i)=reshape(mvnpdf(data,mu(i,:),Sigma{i}),N,1);
        end
 
        for i=1:N
            sum_j=10e-200;
            for j=1:K
                sum_j = sum_j +p(j)*mvns(i,j);
            end
            for j=1:K
                X_grad{j} = X_grad{j} - (p(j)*grad_A(data(i,:)-mu(j,:),X{j},mvns(i,j))/sum_j);
                mu_grad(j,:) = mu_grad(j,:) - p(j)*grad_mu(data(i,:),mu(j,:), X{j}, mvns(i,j))/sum_j;
            end 
            a_grad = a_grad - grad_a_2(p,mvns(i,:))/ sum_j;        
        end 
        grad_struct = struct('X',{X_grad},'mu',mu_grad, 'a', a_grad);     
    end
    
    % Define the problem cost function and its gradient.
    problem.cost  = @negloglikelihood;
    problem.grad = @grad_negloglikelihood;

    %Solve.
    function stats = statsfunc(problem,x,stats)
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
    a0 = exp(w0-max(w0)) ./ sum(exp(w0-max(w0)));
    a0 = reshape(a0,K,1);
    x0 = struct("X",X0,"mu",mu0,"a",a0);
    options=struct("maxiter",maxiter,"minstepsize",1e-12,...
    "memory",10,"statsfun",@statsfunc,"verbosity",1,"tolgradnorm",1e-6,"stopfun",@stopcrit);
    if myrlbfgs
    	[x, xcost, info,options,singular,lsfailed] = my_rlbfgs(problem,x0,options); 
    else
        [x, xcost, info] = rlbfgs(problem,x0,options); 
    end 
    
    costs = [info.cost];
    times = [info.time];
    grads = [info.grad];
    stepsizes = [info.stepsize];
    costevals = [info.costevals];
    
    out=struct("xs",[info.x],"costs",costs,"times",times,"grads",grads,"stepsizes",...
        stepsizes,"gradnorms",[info.gradnorm],"costevals",costevals,...
        "singular",singular,"samecost",samecost,"lsfailed",lsfailed);
end
