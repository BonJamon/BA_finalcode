function [w_grad, S_grad] = grad_negloglikelihood(data,S,w)
    N = size(data,1);
    D = size(data,2);
    K = size(S,1);
    y = [data, ones(N,1)];
    S_grad = cellmat(K,1,D+1,D+1);
    w_grad = zeros(1,K);

    p = exp(w-max(w)) ./ sum(exp(w-max(w)));

    mu=zeros(1,D+1);

    mvns = zeros(N,K);
    for i=1:K
        mvns(:,i)=reshape(mvnpdf(y,mu,S{i}),N,1);
    end

    for i=1:N
        sum_j=0;
        for j=1:K
            sum_j = sum_j +p(j)*mvns(i,j);
        end
        for j=1:K
            S_grad{j} = S_grad{j} - (p(j)*grad_S(y(i,:),S{j},mvns(i,j))/sum_j);
        end 
        for j=1:(K-1)
            w_grad(j) = w_grad(j) - (grad_w(p(j),mvns(i,j)) / sum_j);
        end

    end 
    