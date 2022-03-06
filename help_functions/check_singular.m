function out = check_singular(x)
    if isfield(x,"S")
        Sigma=x.S;
    elseif isfield(x,"A")
        Sigma=get_S_from_A(x.A);
    elseif isfield(x,"X")
        Sigma=get_S_from_A(x.X);
    elseif isfield(x,"Sigma")
        Sigma=x.Sigma;
    end 
    out=false;
    K=size(Sigma,1);
    for i=1:K
        [T,p]=cholcov(Sigma{i},0);
        if p~=0
            out=true;
            break
        end 
    end 
