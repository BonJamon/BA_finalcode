%Script to Compute all results
%Remark: Often copied out: Reparametrized algorithms: Not tested for
%N=25000 or D=20

Ns=[250,2500,25000];
Ks=[2,5];
Ds=[2,5,20];
cs=[0.2,1,5];
seeds=[123,21,34,554,15];
n_sets=size(seeds,2);

maxiter=1500;

mkdir(".", "Results");
mkdir("./Results", "results");

for i=1:size(Ns,2)
    relpathN_store = "./Results/results";
    relpathN_data = "./data";
    Nstr = "N"+int2str(Ns(i));
    mkdir(relpathN_store,Nstr);
    for j=1:size(Ks,2)
        relpathK_store=relpathN_store+"/"+Nstr;
        relpathK_data=relpathN_data+"/"+Nstr;
        Kstr="K"+int2str(Ks(j));
        mkdir(relpathK_store,Kstr);
        for k=1:size(Ds,2)
            relpathD_store=relpathK_store+"/"+Kstr;
            relpathD_data=relpathK_data+"/"+Kstr;
            Dstr="D"+int2str(Ds(k));
            mkdir(relpathD_store,Dstr);
            for l=1:size(cs,2)
                relpathc_store=relpathD_store+"/"+Dstr;
                relpathc_data=relpathD_data+"/"+Dstr;
                if l==1
                    cstr="c"+"02";
                else
                    cstr="c"+int2str(cs(l));
                end 
                mkdir(relpathc_store,cstr);
                for m=1:n_sets
                    if true
                        relpath_store=relpathc_store+"/"+cstr;
                        relpath_data=relpathc_data+"/"+cstr;
                        filename_data=relpath_data+"/"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        data=load(filename_data);
                        data=data.data.data;
                        

                        [w0,mu0,X0]=initialize_Parameters(data,Ks(j),seeds(m));
                        
                        filename_w_man = relpath_store+"/wMan_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if ~isfile(filename_w_man)
                            res_w_man=lbfgs_man_w(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),maxiter,true);
                            save(filename_w_man,"res_w_man");
                        end 
                        
                        filename_w = relpath_store+"/w_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if ~isfile(filename_w)
                            res_w=lbfgs_w(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),maxiter,true);
                            save(filename_w,"res_w");
                        end 
                        
                        filename_a = relpath_store+"/a_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if ~isfile(filename_a)
                            res_a=lbfgs_a(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),maxiter,true);
                            save(filename_a,"res_a");
                        end
                        %{
                        filename_rep = relpath_store+"/rep_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if~isfile(filename_rep)
                            res_rep=lbfgs_reparam(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),false,maxiter,true);
                            save(filename_rep,"res_rep");
                        end 
                        
                        filename_rep_a = relpath_store+"/repA_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if~isfile(filename_rep_a)
                            res_rep_a=lbfgs_reparam_a(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),false,maxiter,true);
                            save(filename_rep_a,"res_rep_a");
                        end 
                        
                        filename_rep_a_pen = relpath_store+"/repAPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if~isfile(filename_rep_a_pen)
                            res_rep_a_pen=lbfgs_reparam_a(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),true,maxiter,true);
                            save(filename_rep_a_pen,"res_rep_a_pen");
                        end 
                        
                        filename_rep_pen = relpath_store+"/repPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if~isfile(filename_rep_pen)
                            res_rep_pen=lbfgs_reparam(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),true,maxiter,true);
                            save(filename_rep_pen,"res_rep_pen");
                        end 
                        %}
                        filename_man = relpath_store+"/man_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if ~isfile(filename_man)
                            res_man=lbfgs_manifold(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),false,maxiter,true);
                            save(filename_man,"res_man");
                        end 
                        
                        filename_man_pen = relpath_store+"/manPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if~isfile(filename_man_pen)
                            res_man_pen=lbfgs_manifold(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),true,maxiter,true);
                            save(filename_man_pen,"res_man_pen");
                        end 
                        
                        filename_EM = relpath_store+"/EM_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if~isfile(filename_EM)
                            res_EM=EM(data,w0,mu0,X0,Ns(i),Ks(j),Ds(k),maxiter);
                            save(filename_EM,"res_EM");
                        end 
                        %{
                        [w0_rand,mu0_rand,X0_rand]=initialize_random(Ks(j),Ds(k),seeds(m));
                        res_rep_penRand=lbfgs_reparam(data,w0_rand,mu0_rand,X0_rand,Ns(i),Ks(j),Ds(k),true,maxiter,true);
                        %res_rep_a_penRand=lbfgs_reparam_a(data,w0_rand,mu0_rand,X0_rand,Ns(i),Ks(j),Ds(k),true,maxiter,true);
                        res_man_penRand=lbfgs_manifold(data,w0_rand,mu0_rand,X0_rand,Ns(i),Ks(j),Ds(k),true,maxiter,true);
                        res_EMRand=EM(data,w0_rand,mu0_rand,X0_rand,Ns(i),Ks(j),Ds(k),maxiter);

                        filename_rep_penRand = relpath_store+"/repPenRand_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        %filename_rep_a_penRand = relpath_store+"/repAPenRand_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        filename_man_penRand = relpath_store+"/manPenRand_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        filename_EMRand = relpath_store+"/EMRand_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        save(filename_rep_penRand,"res_rep_penRand");
                        %save(filename_rep_a_penRand,"res_rep_a_penRand");
                        save(filename_man_penRand,"res_man_penRand");
                        save(filename_EMRand,"res_EMRand");
                        %}
                    end 
                end
            end
        end 
    end
end 