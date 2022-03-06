%Script computing relative amount of times several algorithms produced the
%same cost 
Ns=[250,2500,25000];
Ks=[2,5];
Ds=[2,5,20];
cs=[0.2,1,5];
n_sets=5;

count=0;
for i=1:size(Ns,2)
    relpathN_store = "./Results/results";
    Nstr = "N"+int2str(Ns(i));
    for j=1:size(Ks,2)
        relpathK_store=relpathN_store+"/"+Nstr;
        Kstr="K"+int2str(Ks(j));
        for k=1:size(Ds,2)
            relpathD_store=relpathK_store+"/"+Kstr;
            Dstr="D"+int2str(Ds(k));
            for l=1:size(cs,2)
                relpathc_store=relpathD_store+"/"+Dstr;
                if l==1
                    cstr="c"+"02";
                else
                    cstr="c"+int2str(cs(l));
                end 
                for m=1:n_sets
                    relpath_store=relpathc_store+"/"+cstr;

                    filename_w = relpath_store+"/w_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_w)
                        res_w = load(filename_w).res_w;
                        count=count+check_same_cost(res_w);
                    end 

                    filename_a = relpath_store+"/a_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_a)
                        res_a=load(filename_a).res_a;
                        count=count+check_same_cost(res_a);
                    end

                    filename_rep = relpath_store+"/rep_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_rep)
                        res_rep=load(filename_rep).res_rep;
                        count=count+check_same_cost(res_rep);
                    end 

                    filename_rep_a = relpath_store+"/repA_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_rep_a)
                        res_rep_a=load(filename_rep_a).res_rep_a;
                        count=count+check_same_cost(res_rep_a);
                    end 

                    filename_rep_a_pen = relpath_store+"/repAPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_rep_a_pen)
                        res_rep_a_pen=load(filename_rep_a_pen).res_rep_a_pen;
                        count=count+check_same_cost(res_rep_a_pen);
                    end 

                    filename_rep_pen = relpath_store+"/repPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_rep_pen)
                        res_rep_pen=load(filename_rep_pen).res_rep_pen;
                        count=count+check_same_cost(res_rep_pen);
                    end  

                    filename_man = relpath_store+"/man_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_man)
                        res_man=load(filename_man).res_man;
                        count=count+check_same_cost(res_man);
                    end 

                    filename_man_pen = relpath_store+"/manPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_man_pen)
                        res_man_pen=load(filename_man_pen).res_man_pen;
                        count=count+check_same_cost(res_man_pen);
                    end   

                    filename_EM = relpath_store+"/EM_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                    if isfile(filename_EM)
                        res_EM=load(filename_EM).res_EM;
                        count=count+check_same_cost(res_EM);
                    end  
                end
            end
        end 
    end
end 
disp(count)
function count = check_same_cost(res)
    count=0;
    costs=res.costs;
    if costs(end)==costs(end-1)
        count=1;
    end 
end 