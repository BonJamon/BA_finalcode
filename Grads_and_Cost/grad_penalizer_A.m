function grad = grad_penalizer_A(A,Psi)
    p=0.01;
    b=1;
    invS = inv(A*A');
    grad=-invS*A*p+0.5*b*invS*(Psi+Psi')*invS*A;