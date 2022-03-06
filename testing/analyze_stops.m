function fig = analyze_stops()
    %Counts why algorithms stopped

    check_reparam=false;
    %settings
    tolgradnorm=1e-6;
    minstepsize=1e-12;
    maxiter=1500;
    
    %init
    if check_reparam
        stats=zeros(10,7);
        [k_w,k_a,k_rep,k_rep_a,k_rep_pen,k_rep_a_pen,k_man,k_man_pen,k_EM,k_w_man]=deal(0);
    else 
        stats=zeros(6,7);
        [k_w,k_a,k_man,k_man_pen,k_EM,k_w_man]=deal(0);
    end 
    Ns=[250,2500,25000];
    Ks=[2,5];
    Ds=[2,5,20];
    cs=[0.2];
    n_sets=5;
    
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
                        if check_reparam
                            statscurrent=zeros(10,7);
                        else 
                            statscurrent=zeros(6,7);
                        end
                        filename_w = relpath_store+"/w_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_w)
                            k_w=k_w+1;
                            res_w = load(filename_w).res_w;
                            res_w=adapt_res(res_w,tolgradnorm,minstepsize,maxiter);
                            statscurrent(1,:)=[res_w.gradnorms(end),res_w.gradnorms(end),res_w.stepsizes(end),...
                                res_w.singular,res_w.samecost,isnan(res_w.costs(end)),res_w.lsfailed];
                        end 
                      
                        filename_a = relpath_store+"/a_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_a)
                            k_a=k_a+1;
                            res_a=load(filename_a).res_a;
                            res_a=adapt_res(res_a,tolgradnorm,minstepsize,maxiter);
                            statscurrent(2,:)=[res_a.gradnorms(end),res_a.gradnorms(end),res_a.stepsizes(end),...
                                res_a.singular,res_a.samecost,isnan(res_a.costs(end)),res_a.lsfailed];
                        end
                        if check_reparam
                            filename_rep = relpath_store+"/rep_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep)
                                k_rep=k_rep+1;
                                res_rep=load(filename_rep).res_rep;
                                res_rep=adapt_res(res_rep,tolgradnorm,minstepsize,maxiter);
                                statscurrent(3,:)=[res_rep.gradnorms(end),res_rep.gradnorms(end),res_rep.stepsizes(end),...
                                    res_rep.singular,res_rep.samecost,isnan(res_rep.costs(end)),res_rep.lsfailed];
                            end 

                            filename_rep_a = relpath_store+"/repA_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep_a)
                                k_rep_a=k_rep_a+1;
                                res_rep_a=load(filename_rep_a).res_rep_a;
                                res_rep_a=adapt_res(res_rep_a,tolgradnorm,minstepsize,maxiter);
                                statscurrent(4,:)=[res_rep_a.gradnorms(end),res_rep_a.gradnorms(end),res_rep_a.stepsizes(end),...
                                    res_rep_a.singular,res_rep_a.samecost,isnan(res_rep_a.costs(end)),res_rep_a.lsfailed];
                            end 

                            filename_rep_a_pen = relpath_store+"/repAPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep_a_pen)
                                k_rep_a_pen=k_rep_a_pen+1;
                                res_rep_a_pen=load(filename_rep_a_pen).res_rep_a_pen;
                                res_rep_a_pen=adapt_res(res_rep_a_pen,tolgradnorm,minstepsize,maxiter);
                                statscurrent(5,:)=[res_rep_a_pen.gradnorms(end),res_rep_a_pen.gradnorms(end),res_rep_a_pen.stepsizes(end),...
                                    res_rep_a_pen.singular,res_rep_a_pen.samecost,isnan(res_rep_a_pen.costs(end)),res_rep_a_pen.lsfailed];
                            end 

                            filename_rep_pen = relpath_store+"/repPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep_pen)
                                k_rep_pen=k_rep_pen+1;
                                res_rep_pen=load(filename_rep_pen).res_rep_pen;
                                res_rep_pen=adapt_res(res_rep_pen,tolgradnorm,minstepsize,maxiter);
                                statscurrent(6,:)=[res_rep_pen.gradnorms(end),res_rep_pen.gradnorms(end),res_rep_pen.stepsizes(end),...
                                    res_rep_pen.singular,res_rep_pen.samecost,isnan(res_rep_pen.costs(end)),res_rep_pen.lsfailed];
                            end  
                         
                        
                            filename_man = relpath_store+"/man_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_man)
                                k_man=k_man+1;
                                res_man=load(filename_man).res_man;
                                res_man=adapt_res(res_man,tolgradnorm,minstepsize,maxiter);
                                statscurrent(7,:)=[res_man.gradnorms(end),res_man.gradnorms(end),res_man.stepsizes(end),...
                                    res_man.singular,res_man.samecost,isnan(res_man.costs(end)),res_man.lsfailed];
                            end 

                            filename_man_pen = relpath_store+"/manPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_man_pen)
                                k_man_pen=k_man_pen+1;
                                res_man_pen=load(filename_man_pen).res_man_pen;
                                res_man_pen=adapt_res(res_man_pen,tolgradnorm,minstepsize,maxiter);
                                statscurrent(8,:)=[res_man_pen.gradnorms(end),res_man_pen.gradnorms(end),res_man_pen.stepsizes(end),...
                                    res_man_pen.singular,res_man_pen.samecost,isnan(res_man_pen.costs(end)),res_man_pen.lsfailed];
                            end   

                            filename_EM = relpath_store+"/EM_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_EM)
                                k_EM=k_EM+1;
                                disp(filename_EM)
                                res_EM=load(filename_EM).res_EM;
                                res_EM=adapt_resEM(res_EM,maxiter);
                                statscurrent(9,:)=[0,1,1,res_EM.singular,0,isnan(res_EM.costs(end)),0];
                            end 
                            filename_w_man = relpath_store+"/wMan_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_w_man)
                                k_w_man=k_w_man+1;
                                res_wMan = load(filename_w_man).res_w_man;
                                res_wMan=adapt_res(res_wMan,tolgradnorm,minstepsize,maxiter);
                                statscurrent(10,:)=[res_wMan.gradnorms(end),res_wMan.gradnorms(end),res_wMan.stepsizes(end),...
                                    res_wMan.singular,res_wMan.samecost,isnan(res_wMan.costs(end)),res_wMan.lsfailed];
                            end 
                        else 
                            filename_man = relpath_store+"/man_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_man)
                                k_man=k_man+1;
                                res_man=load(filename_man).res_man;
                                res_man=adapt_res(res_man,tolgradnorm,minstepsize,maxiter);
                                statscurrent(3,:)=[res_man.gradnorms(end),res_man.gradnorms(end),res_man.stepsizes(end),...
                                    res_man.singular,res_man.samecost,isnan(res_man.costs(end)),res_man.lsfailed];
                            end 

                            filename_man_pen = relpath_store+"/manPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_man_pen)
                                k_man_pen=k_man_pen+1;
                                res_man_pen=load(filename_man_pen).res_man_pen;
                                res_man_pen=adapt_res(res_man_pen,tolgradnorm,minstepsize,maxiter);
                                statscurrent(4,:)=[res_man_pen.gradnorms(end),res_man_pen.gradnorms(end),res_man_pen.stepsizes(end),...
                                    res_man_pen.singular,res_man_pen.samecost,isnan(res_man_pen.costs(end)),res_man_pen.lsfailed];
                            end   

                            filename_EM = relpath_store+"/EM_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_EM)
                                k_EM=k_EM+1;
                                disp(filename_EM)
                                res_EM=load(filename_EM).res_EM;
                                res_EM=adapt_resEM(res_EM,maxiter);
                                statscurrent(5,:)=[0,1,1,res_EM.singular,0,isnan(res_EM.costs(end)),0];
                            end 
                            filename_w_man = relpath_store+"/wMan_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_w_man)
                                k_w_man=k_w_man+1;
                                res_wMan = load(filename_w_man).res_w_man;
                                res_wMan=adapt_res(res_wMan,tolgradnorm,minstepsize,maxiter);
                                statscurrent(6,:)=[res_wMan.gradnorms(end),res_wMan.gradnorms(end),res_wMan.stepsizes(end),...
                                    res_wMan.singular,res_wMan.samecost,isnan(res_wMan.costs(end)),res_wMan.lsfailed];
                            end 
                        end
                        if isfile(filename_EM)
                            statscurrent(:,2)=statscurrent(:,2)<tolgradnorm;
                            statscurrent(:,3)=statscurrent(:,3)<minstepsize;
                            stats=stats+statscurrent;
                        end 

                    end
                end
            end 
        end
    end 
    if check_reparam
        ks=[k_w,k_a,k_rep,k_rep_a,k_rep_a_pen,k_rep_pen,k_man,k_man_pen,k_EM,k_w_man];
    else 
        ks=[k_w,k_a,k_man,k_man_pen,k_EM,k_w_man];
    end
    ks=reshape(ks,[],1);
    stats=stats./ks;
    stats(:,2:7)=100*stats(:,2:7);
    
    
    varNames=["Mean gradnorm","% tolgradnorm","% minstepsize",...
        "% singular","% same cost","% ~posdef","% LS failed"];
    
    out_table=array2table(stats,"VariableNames",varNames);
    if check_reparam
        Algorithms=["w","a","rep","repA","repAPen","repPen","man","manPen","EM","wMan"];
    else 
        Algorithms=["w","a","man","manPen","EM","wMan"];
    end 
    Algorithms=Algorithms.reshape([],1);
    out_table=[table(Algorithms) out_table];
    
    fig = uifigure("Position",[120,120,760,360]);
    out_table=uitable(fig,"Data",out_table,"Position",[20,20,720,320]);
    
end 
function res = adapt_res(res,tolgradnorm,minstepsize,maxiter)
    %Implement: Should find when tolgradnorm would have triggered even when
    %minstepsize is
    
    %Find index where it would have stopped
    t=find(res.gradnorms<tolgradnorm);
    if isempty(t)
        t=maxiter;
    else 
        t=t(1);
    end 
    %{
    s=find(res.stepsizes<minstepsize);
    if isempty(s)
        s=maxiter;
    else 
        s=s(1);
    end 
    %}
    iters=numel(res.costs);
    k=min([t+1,maxiter+1,iters]);
    %Change the structure
    names=fieldnames(res);
    s=size(names,1);
    names=names(1:s-3);
    if k<iters
        res.singular=false;
        res.lsfailed=false;
        res.samecost=false;
    end
    for i=1:numel(names)
        name_i=string(names(i));
        res.(name_i)=res.(name_i)(1:k);
    end 
end 
function res=adapt_resEM(res,maxiter)
    names=fieldnames(res);
    iters=numel(res.costs);
    if maxiter+1<iters
        res.singular=false;
    end 
    k=min([maxiter,iters]);
    s=size(names,1);
    names=names(1:s-1);
    for i=1:numel(names)
        name_i=string(names(i));
        if ~isempty(res.(name_i))
            res.(name_i)=res.(name_i)(1:k);
        end 
    end 
end
    