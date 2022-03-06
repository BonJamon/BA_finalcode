%N:Number of Points
%K:Number of Classes
%D:Number of Dimensions
%c:Degree of seperation
Ns=[250,2500,25000];
Ks=[2,5];
Ds=[2,5,20];
cs=[0.2,1,5];
seeds=[1,2,3,4,5];
n_sets=size(seeds,2);

mkdir(".", "data");

for i=1:size(Ns,2)
    relpathN = "./data";
    Nstr = "N"+int2str(Ns(i));
    mkdir(relpathN,Nstr);
    for j=1:size(Ks,2)
        relpathK=relpathN+"/"+Nstr;
        Kstr="K"+int2str(Ks(j));
        mkdir(relpathK,Kstr);
        for k=1:size(Ds,2)
            relpathD=relpathK+"/"+Kstr;
            Dstr="D"+int2str(Ds(k));
            mkdir(relpathD,Dstr);
            for l=1:size(cs,2)
                relpathc=relpathD+"/"+Dstr;
                if l==1
                    cstr="c"+"02";
                else
                    cstr="c"+int2str(cs(l));
                end 
                mkdir(relpathc,cstr);
                for m=1:n_sets
                    relpath=relpathc+"/"+cstr;
                    data = create_centers(Ns(i),Ks(j),Ds(k),cs(l),seeds(m));
                    filename=relpath+"/"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    save(filename,"data");
                end
            end
        end 
    end
end 