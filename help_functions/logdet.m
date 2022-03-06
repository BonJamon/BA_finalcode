function res = logdet(S)
   %S is KxDxD
   %res is Kx1
   K = size(S,1);
   res = zeros(K);
   for i = 1:K
       %Problem: Eigenvalues sometimes 0
       U = chol(S{i});
       res(i)=2*sum(log(diag(U)));
   end
   