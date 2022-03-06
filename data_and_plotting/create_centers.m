function out = create_centers(N,K,D,c,seed)
    rng(seed)
    covs=cellmat(K,1,D,D);
    traces=zeros(K,1);
    i=1;
    while i~=K+1
        X = randn(D,D);
        cov=0.2*(X*X');
        if ~any(eig(cov)<1e-4)
            covs{i}=cov;
            traces(i)=trace(covs{i});
            i=i+1;
        end 
    end 

    % create c-separated cluster centers
    trials = 0;

    x = zeros(D, K);
    x(end, :) = 1;
    while 1 
      rd = randn(D, K);
      rd(end, :) = 0;
      nd = sqrt(sum(rd.^2, 1));
      rd = bsxfun(@rdivide, rd, nd);

      t = sqrt(c) * sqrt(K) * trials/10 * rand(1,K);
      %Problem: Muss ensuren, dass ein bestimmter Abstand der einzelnen
      %Cluser möglich ist. 
      M = D*D*K*(bsxfun(@times, x, cos(t)) + bsxfun(@times, rd, sin(t)));
      
      % check degree of separation
      error = 0;
      for i = 1:K-1
        for j = i+1:K
          if norm(M(:,i)-M(:,j)) < c * max(traces(i),traces(j))
            error = 1;
          end
        end
      end
      if ~error
        break;
      end
      trials = trials + 0.1;
    end
    M=M';
    data=zeros(N,D);
    NpC = floor(N / K);
    N_done=0;
    for i = 1:K-1
        data(N_done+1:N_done+NpC,:)=mvnrnd(M(i,:),covs{i},NpC);
        N_done=N_done+NpC;
    end 
    data(N_done+1:N,:)=mvnrnd(M(K,:),covs{K},N-N_done);
    out=struct("data",data,"Sigma",covs,"mu",M);
    %plot_data(data,M,covs,1/K*ones(1,K));