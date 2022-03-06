function res = log_q(y,S,logdetS)
%y is NxD+1, S is KxD+1xD+1, res is NxK
    N = size(y,1);
    K = size(S,1);
    res = zeros(N,K);

    for i = 1:N
        for j = 1:K
            res(i,j)=logdetS(j) + y(i,:)/S{j}*y(i,:)';
        end
    end
    res = (-0.5)*res;