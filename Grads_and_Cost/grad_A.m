function grad = grad_A(datapoint,A,mvn)
    D=size(datapoint);
    datapoint = reshape(datapoint, [D,1]);
    S=A*A';
    %Z=inv(S);
    grad=mvn*(S\(datapoint'*datapoint)/S*A-S\A);



