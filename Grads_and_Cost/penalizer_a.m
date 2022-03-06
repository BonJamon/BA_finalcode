function out = penalizer_a(a)
    %set c=1 --> just product of array elements
    K=size(a,2);
    out=reshape(a,K,1);
    out=prod(out);