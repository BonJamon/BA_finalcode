function res = logsumexp(alpha, log_q)
%alpha:K log_q: NxK
    N = size(log_q, 1);
    K=size(alpha);
    res = zeros(N,1);
    %c=max(log_q());-->need logsumexptrick
    for i =1:N
        for j=1:K
            log_q(i,:)= exp(log_q(i,:)+log(alpha));
        end
    res=log(sum(log_q,2));
    end
    