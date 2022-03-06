function res = get_A_from_S(S)
    K=size(S,1);
    res = cell(K,1);
    for i=1:K
        res{i,:}=chol(S{i},"lower");
    end