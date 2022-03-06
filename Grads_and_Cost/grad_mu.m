function grad = grad_mu(datapoint,mu,X,mvn)
    %gradient wRt mu when given X (Sigma=XX')
    dmu=datapoint-mu;
    grad = mvn*(dmu/(X*X'));