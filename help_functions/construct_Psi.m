function Psi = construct_Psi(a,b,Delta,k,lambd)
    D=size(lambd,2);
    Psi = zeros(D+1,D+1);
    Psi(1:D,1:D)=a/b*Delta+k*(lambd'*lambd);
    Psi(1:D,D+1)=k*lambd;
    Psi(D+1,1:D)=k*lambd';
    Psi(D+1,D+1)=k;
    

    
