function gradient_checking()
    %General script to test parts or full gradients of functions using
    %checkgradient() of manopt
    N=200;
    K=2;
    D=2;

    c=1;
    data=generate_data(N,c,142);
    datapoint=data(1,:);
    
    %CHECK: Gradient Normal Distribution w.R.t. A
    problem1.M=euclideanfactory(D+1,D+1);
    function res = pN_A(A)
        y=[datapoint,1];
        mu=zeros(1,D+1);
        Sigma=A*A';
        res=mvnpdf(y,mu,Sigma);
    end 
    function grad=grad_pN_A(A)
        y=[datapoint,1];
        mu=zeros(1,D+1);
        Sigma=A*A';
        mvn=mvnpdf(y,mu,Sigma);
        grad=grad_A(y,A,mvn);
    end
    problem1.cost=@pN_A;
    problem1.grad=@grad_pN_A;
    disp("check grad A");
    figure;
    checkgradient(problem1);
    
    %CHECK: w-Term w.R.t. w
    problem2.M=euclideanfactory(1,K-1);
    mvns_dummy=rand(N,K);
    terms=mvns_dummy(1,:);
    function res = term_w(w)
        w=[w,0];
        alpha=exp(w-max(w))./sum(exp(w-max(w)));
        res=sum(alpha.*terms);
    end 
    function grad = grad_term_w(w)
        w=[w,0];
        alpha=exp(w-max(w))./sum(exp(w-max(w)));
        grad=zeros(1,K-1);
        for i=1:K-1
            grad(i)=alpha(i)*(terms(i)-sum(alpha.*terms));
        end 
    end 
    problem2.cost=@term_w;
    problem2.grad=@grad_term_w;
    disp("check grad w");
    figure;
    checkgradient(problem2);
            
    
    %Check everything using reparam_w, but without penalizer
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
        res=res-penalizer_w(w);
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
            sum_j=0;
            for j=1:K
                sum_j = sum_j + p(j)*mvns(i,j);
            end
            for j=1:K
                A_grad{j} = A_grad{j} - (p(j)*grad_A(y(i,:),A{j},mvns(i,j))/sum_j);
            end 
            for j=1:(K-1)
                w_grad(j) = w_grad(j) - p(j)*(mvns(i,j)-sum(p.*mvns(i,:)))/ sum_j;
            end
        end 
        w_grad = w_grad - grad_penalizer_w(w);
        grad_struct = struct('A',{A_grad}, 'w', w_grad);     
    end
    problem.cost  = @negloglikelihood;
    problem.grad = @grad_negloglikelihood;
    disp("check gradient w reparam");
    figure;
    checkgradient(problem)
    
    %CHECK: penalizer_w
    problem4.M=euclideanfactory(1,K-1);
    problem4.cost=@penalizer_w;
    problem4.grad=@grad_penalizer_w;
    disp("check gradient penalizer_w");
    figure;
    checkgradient(problem4);
    
    %TODO!
    %CHECK: penalizer_S
    a=1;
    b=1;
    Delta=0.01*cov(data);
    k=1;
    lambd=mean(data);
    Psi=construct_Psi(a,b,Delta,k,lambd);
    problem5.M=sympositivedefinitefactory(D+1);
    function res = pen_S(S)
        res=penalizer_S(S,Psi);
    end 
    function res =grad_pen_S(S)
        res=grad_penalizer_S(S,Psi);
    end 
    problem5.cost=@pen_S;
    problem5.grad=@grad_pen_S;
    disp("check penalizer S");
    figure;
    checkgradient(problem5);
    
    %CHECK: Gradient reparam_a
    struct1 = struct('A', powermanifold(euclideanfactory(D+1,D+1),K), 'a', multinomialfactory(K,1));
    manifold = productmanifold(struct1);
    problem6.M = manifold;

    function res = negloglikelihood_a(varstruct)  
        A = varstruct.A;
        alpha = reshape(varstruct.a,1,K);
        S=get_S_from_A(A);
        
        N = size(data,1);
        y = [data, ones(N,1)];
        
        mu=zeros(K,D+1);
        res=negloglikelihood_Sigma(y,alpha,mu,S,true);
    end

    function grad_struct = grad_negloglikelihood_a(varstruct)
        A = varstruct.A;
        p=varstruct.a;
        y = [data, ones(N,1)];
  
        A_grad = cellmat(K,1,D+1,D+1);
        a_grad = zeros(K,1);
     

        mu=zeros(1,D+1);

        mvns = zeros(N,K);
        for i=1:K
            mvns(:,i)=reshape(mvnpdf(y,mu,A{i}*A{i}'),N,1);
        end
 
        for i=1:N
            sum_j=10e-200;
            for j=1:K
                sum_j = sum_j +p(j)*mvns(i,j);
            end
            for j=1:K
                A_grad{j} = A_grad{j} - (p(j)*grad_A(y(i,:),A{j},mvns(i,j))/sum_j);
            end 
            a_grad=a_grad-grad_a_2(p,mvns(i,:))/sum_j;
                
        end
        grad_struct = struct('A',{A_grad}, 'a', a_grad);     
    end
    
    % Define the problem cost function and its gradient.
    problem6.cost  = @negloglikelihood_a;
    problem6.grad = @grad_negloglikelihood_a;
    disp("check reparam_a")
    figure;
    checkgradient(problem6);
    
    problem7.M=multinomialfactory(K,1);
    terms=mvns_dummy(1,:);
    %CHECK: Gradient a
    function res=inner_a(a)
        a=reshape(a,1,K);
        res=sum(terms.*a);
    end 
    function grad_riem = grad_inner_a(a)
        grad_euclid=reshape(terms,K,1);
        Z=grad_euclid.*a;
        alpha=ones(1,K)*Z;
        grad_riem=Z-(ones(K,1)*alpha).*a;
    end 
    problem7.cost=@inner_a;
    problem7.grad=@grad_inner_a;
    %problem7=manoptAD(problem7);
    disp("check inner a");
    figure;
    checkgradient(problem7);
    
    %CHECK: Gradient Normal Distribution w.R.t. S
    problem8.M=sympositivedefinitefactory(D+1);
    function res = pN_S(S)
        y=[datapoint,1];
        mu=zeros(1,D+1);
        res=mvnpdf(y,mu,S);
    end 
    
    function grad=grad_pN_S(S)
        y=[datapoint,1];
        mu=zeros(1,D+1);
        mvn=mvnpdf(y,mu,S);
        grad=grad_S(y,S,mvn);
    end

    problem8.cost=@pN_S;
    problem8.grad=@grad_pN_S;
    disp("check grad S");
    figure;
    checkgradient(problem8);
    
    %CHECK: Penalizer_A
    problem9.M=euclideanfactory(D+1,D+1);
    function x=pen_A(A)
        x=penalizer_S(A*A',Psi);
    end 
    function x=grad_pen_A(A)
        x=grad_penalizer_A(A,Psi);
    end 
    problem9.cost=@pen_A;
    problem9.grad=@grad_pen_A;
    
    disp("check penalizer A");
    figure;
    checkgradient(problem9);
    X4=randn(D,D);
    %CHECK: Grad Normal Distribution w.R.t. mu
    problem10.M=euclideanfactory(1,D);
    function res = cost_mu(mu)
        res = mvnpdf(data(1,:),mu,X4*X4');
    end 
    function res = grad_mu1(mu)
        mvn=mvnpdf(data(1,:),mu,X4*X4');
        res=grad_mu(data(1,:),mu,X4,mvn);
    end 
    problem10.cost=@cost_mu;
    problem10.grad=@grad_mu1;
    figure;
    checkgradient(problem10);
end
