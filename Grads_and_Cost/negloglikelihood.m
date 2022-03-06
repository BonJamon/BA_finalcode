function res = negloglikelihood(data, S, w)
        N = size(data,1);
        logdetS = logdet(S);

        y = [data, ones(N,1)];

        alpha= exp(w);
        alpha = alpha ./ sum(alpha);
        
        loglikevec = logsumexp(alpha,log_q(y,S,logdetS));

        res = -sum(loglikevec);
    
    