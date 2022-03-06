function out = penalizer_S(S,Psi)
    p=0.01;
    b=1;
    logdetS=log(det(S));
    out = (-0.5)*p*logdetS-0.5*b*trace(Psi*inv(S));
    %out=(-0.5)*b*trace(Psi/S);
    
