function grad = grad_w(p,mvn)
%mvn: density value, p:value of corresponding alpha
%INCORRECT!
    grad=p.*(mvn-sum(p.*mvn));