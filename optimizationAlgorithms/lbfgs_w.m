function out = lbfgs_w(data,w0,mu0,X0,N,K,D,maxiter,myrlbfgs)
    %Create Problem Structure
    struct1 = struct('X', powermanifold(euclideanfactory(D,D),K),'mu',euclideanfactory(K,D), 'w', euclideanfactory(1,K-1));
    manifold = productmanifold(struct1);
    problem.M = manifold;
    singular=false;
    samecost=false;
    lsfailed=false;

    function res = negloglikelihood(varstruct) 
        X = varstruct.X;
        mu = varstruct.mu;
        w = varstruct.w;
        Sigma = get_S_from_A(X);
        res = negloglikelihood_Sigma(data,w,mu,Sigma,false);
    end

    function grad_struct = grad_negloglikelihood(varstruct)
        X = varstruct.X;
        mu = varstruct.mu;
        w = varstruct.w;
        w = [w,0];
        
        X_grad = cellmat(K,1,D,D);
        mu_grad = zeros(K,D);
        w_grad = zeros(1,K-1);
        
        Sigma=get_S_from_A(X);
        p = exp(w-max(w)) ./ sum(exp(w-max(w)));

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
            for j=1:(K-1)
                w_grad(j) = w_grad(j) - p(j)*(mvns(i,j)-sum(p.*mvns(i,:))) / sum_j;
            end
          
                
        end 
        grad_struct = struct('X',{X_grad},'mu',mu_grad, 'w', w_grad);     
    end
    
    % Define the problem cost function and its gradient.
    problem.cost  = @negloglikelihood;
    problem.grad = @grad_negloglikelihood;

    %Solve.
    w0 = w0(1:(K-1));
    x0 = struct("X",X0,"mu",mu0,"w",w0);
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
    options=struct("maxiter",maxiter,"minstepsize",1e-12,"memory",10,...
        "statsfun",@statsfunc,"verbosity",1,"tolgradnorm",1e-6,"stopfun",@stopcrit);
    if myrlbfgs
    	[x, xcost, info,options,singular,lsfailed] = my_rlbfgs(problem,x0,options); 
    else
        [x, xcost, info] = rlbfgs(problem,x0,options); 
    end 
    grads=[info.grad];
    stepsizes=[info.stepsize];
    costevals=[info.costevals];
    out = struct("costs", [info.cost],"times",[info.time],"grads",grads,...
        "stepsizes",stepsizes,"gradnorms",[info.gradnorm],"costevals",costevals,...
        "singular",singular,"samecost",samecost,"lsfailed",lsfailed);

    
end
