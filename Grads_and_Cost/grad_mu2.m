function grad = grad_mu2(datapoint,mu,Sigma,mvn)
%gradient wRt mu when given Sigma
    dmu=datapoint-mu;
    grad = mvn*(dmu/Sigma);