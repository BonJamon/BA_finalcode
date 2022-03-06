function grad = grad_penalizer_S(S,Psi)
    p=0.01;
    b=1;
    grad=(-0.5)*p*S+0.5*b*Psi;