function res = get_S_from_A(A)
    K=size(A,1);
    res = cell(K,1);
    for i=1:K
        res{i,:}=A{i}*A{i}';
    end