function grad = grad_S(datapoint,S,mvn)
    D=size(datapoint);
    datapoint = reshape(datapoint, [D,1]);
    Z=inv(S);
    grad=0.5*mvn*(datapoint'*datapoint-S);