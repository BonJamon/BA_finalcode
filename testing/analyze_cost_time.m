function fig = analyze_cost_time()
    %Computes relative or absolute cost, time, iteration values over all
    %results
    %Just a general script: The implementation of what actually is computed
    %varies
    check_reparam=false;

    %settings
    tolgradnorm=1e-6;
    minstepsize=1e-12;
    maxiter=1500;
    
    %init
    count=0;
    Ns=[250,2500,25000];
    Ks=[2,5];
    Ds=[2,5,20];
    cs=[0.2,1,5];
    n_sets=5;
    
    n_results=size(Ns,2)*size(Ks,2)*size(Ds,2)*size(cs,2)*n_sets;
    costs_w=zeros(n_results,1);
    times_w=zeros(n_results,1);
    costs_a=zeros(n_results,1);
    times_a=zeros(n_results,1);
    if check_reparam
        costs_rep=zeros(n_results,1);
        times_rep=zeros(n_results,1);
        costs_repA=zeros(n_results,1);
        times_repA=zeros(n_results,1);
        costs_repAPen=zeros(n_results,1);
        times_repAPen=zeros(n_results,1);
        costs_repPen=zeros(n_results,1);
        times_repPen=zeros(n_results,1);
    end 
    costs_man=zeros(n_results,1);
    times_man=zeros(n_results,1);
    costs_manPen=zeros(n_results,1);
    times_manPen=zeros(n_results,1);
    costs_EM=zeros(n_results,1);
    times_EM=zeros(n_results,1);
    costs_wMan=zeros(n_results,1);
    times_wMan=zeros(n_results,1);
    
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
                    if cs(l)==0.2
                        cstr="c"+"02";
                    else
                        cstr="c"+int2str(cs(l));
                    end 
                    for m=1:n_sets
                        relpath_store=relpathc_store+"/"+cstr;
                        filename_w = relpath_store+"/w_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_w)
                            res_w = load(filename_w).res_w;
                            res_w=adapt_res(res_w,tolgradnorm,minstepsize,maxiter);
                        end 
                      
                        filename_a = relpath_store+"/a_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_a)
                            res_a=load(filename_a).res_a;
                            res_a=adapt_res(res_a,tolgradnorm,minstepsize,maxiter);
                        end
                        
                        if check_reparam
                            filename_rep = relpath_store+"/rep_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep)
                                res_rep=load(filename_rep).res_rep;
                                res_rep=adapt_res(res_rep,tolgradnorm,minstepsize,maxiter);
                            end 

                            filename_rep_a = relpath_store+"/repA_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep_a)
                                res_rep_a=load(filename_rep_a).res_rep_a;
                                res_rep_a=adapt_res(res_rep_a,tolgradnorm,minstepsize,maxiter);
                            end 

                            filename_rep_a_pen = relpath_store+"/repAPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep_a_pen)
                                res_rep_a_pen=load(filename_rep_a_pen).res_rep_a_pen;
                                res_rep_a_pen=adapt_res(res_rep_a_pen,tolgradnorm,minstepsize,maxiter);
                            end 

                            filename_rep_pen = relpath_store+"/repPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                            if isfile(filename_rep_pen)
                                res_rep_pen=load(filename_rep_pen).res_rep_pen;
                                res_rep_pen=adapt_res(res_rep_pen,tolgradnorm,minstepsize,maxiter);
                            end  
                        end
                        filename_man = relpath_store+"/man_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_man)
                            res_man=load(filename_man).res_man;
                            res_man=adapt_res(res_man,tolgradnorm,minstepsize,maxiter);
                        end 
                        
                        filename_man_pen = relpath_store+"/manPen_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_man_pen)
                            res_man_pen=load(filename_man_pen).res_man_pen;
                            res_man_pen=adapt_res(res_man_pen,tolgradnorm,minstepsize,maxiter);
                        end   
                        
                        filename_EM = relpath_store+"/EM_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_EM)
                            count=count+1;
                            disp(filename_EM)
                            res_EM=load(filename_EM).res_EM;
                            res_EM=adapt_resEM(res_EM,maxiter);
                        end 
                        
                        filename_wMan = relpath_store+"/wMan_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_wMan)
                            res_wMan = load(filename_wMan).res_w_man;
                            res_wMan=adapt_res(res_wMan,tolgradnorm,minstepsize,maxiter);
                        end 
                        
                        costs_w(count)=1;
                        times_w(count)=1;
                        costs_a(count)=res_a.costs(end)/res_w.costs(end);
                        times_a(count)=res_a.times(end)/res_w.times(end);
                        if check_reparam
                            costs_rep(count)=res_rep.costs(end)/res_w.costs(end);
                            times_rep(count)=res_rep.times(end)/res_w.times(end);
                            costs_repA(count)=res_rep_a.costs(end)/res_w.costs(end);
                            times_repA(count)=res_rep_a.times(end)/res_w.times(end);
                            costs_repAPen(count)=res_rep_a_pen.costs(end)/res_w.costs(end);
                            times_repAPen(count)=res_rep_a_pen.times(end)/res_w.times(end);
                            costs_repPen(count)=res_rep_pen.costs(end)/res_w.costs(end);
                            times_repPen(count)=res_rep_pen.times(end)/res_w.times(end);
                        end
                        costs_man(count)=1;%res_man.costs(end)/res_w.costs(end);
                        times_man(count)=1;%res_man.times(end)/res_w.times(end);
                        costs_manPen(count)=res_man_pen.costs(end)/res_man.costs(end);%res_w.costs(end);
                        times_manPen(count)=res_man_pen.times(end)/res_man.times(end);%res_w.times(end);
                        costs_EM(count)=res_EM.costs(end)/res_w.costs(end);
                        %been a mistake in EM code for singular case
                        if isempty(res_EM.times)
                            times_EM(count)=NaN;
                        else 
                            times_EM(count)=res_EM.times(end)/res_w.times(end);
                        end 
                        costs_wMan(count)=res_wMan.costs(end)/res_w.costs(end);
                        times_wMan(count)=res_wMan.times(end)/res_w.times(end);
                    end
                end
            end 
        end
    end 
    %For the case not all results are there
    costs_w=costs_w(1:count);
    times_w=times_w(1:count);
    costs_a=costs_a(1:count);
    times_a=times_a(1:count);
    if check_reparam
        costs_rep=costs_rep(1:count);
        times_rep=times_rep(1:count);
        costs_repA=costs_repA(1:count);
        times_repA=times_repA(1:count);
        costs_repAPen=costs_repAPen(1:count);
        times_repAPen=times_repAPen(1:count);
        costs_repPen=costs_repPen(1:count);
        times_repPen=times_repPen(1:count);
    end 
    
    costs_man=costs_man(1:count);
    times_man=times_man(1:count);
    costs_manPen=costs_manPen(1:count);
    times_manPen=times_manPen(1:count);
    costs_EM=costs_EM(1:count);
    times_EM=times_EM(1:count);
    costs_wMan=costs_wMan(1:count);
    times_wMan=times_wMan(1:count);
    
    %Only take costs not NaN
    times_EM=times_EM(~isnan(times_EM));
    if check_reparam
        costs_rep=costs_rep(~isnan(costs_rep));
        costs_repPen=costs_repPen(~isnan(costs_repPen));
        costs_repA=costs_repA(~isnan(costs_repA));
        costs_repAPen=costs_repAPen(~isnan(costs_repAPen));
    end
    costs_man=costs_man(~isnan(costs_man));
    costs_manPen=costs_manPen(~isnan(costs_manPen));
    costs_wMan=costs_wMan(~isnan(costs_wMan));
    
    stats(1,:)=[mean(costs_w),std(costs_w),mean(times_w),std(times_w)];
    stats(2,:)=[mean(costs_a),std(costs_a),mean(times_a),std(times_a)];
    if check_reparam
        stats(3,:)=[mean(costs_rep),std(costs_rep),mean(times_rep),std(times_rep)];
        stats(4,:)=[mean(costs_repA),std(costs_repA),mean(times_repA),std(times_repA)];
        stats(5,:)=[mean(costs_repAPen),std(costs_repAPen),mean(times_repAPen),std(times_repAPen)];
        stats(6,:)=[mean(costs_repPen),std(costs_repPen),mean(times_repPen),std(times_repPen)];
        stats(7,:)=[mean(costs_man),std(costs_man),mean(times_man),std(times_man)];
        stats(8,:)=[mean(costs_manPen),std(costs_manPen),mean(times_manPen),std(times_manPen)];
        stats(9,:)=[mean(costs_EM),std(costs_EM),mean(times_EM),std(times_EM)];
        stats(10,:)=[mean(costs_wMan),std(costs_wMan),mean(times_wMan),std(times_wMan)];
    else 
        stats(3,:)=[mean(costs_man),std(costs_man),mean(times_man),std(times_man)];
        stats(4,:)=[mean(costs_manPen),std(costs_manPen),mean(times_manPen),std(times_manPen)];
        stats(5,:)=[mean(costs_EM),std(costs_EM),mean(times_EM),std(times_EM)];
        stats(6,:)=[mean(costs_wMan),std(costs_wMan),mean(times_wMan),std(times_wMan)];
    end 
        
    
    
    varNames=["Ratio to w cost","std of Ratio to w cost","Mean Ratio to w time","std of Ratio to w time"];
    
    out_table=array2table(stats,"VariableNames",varNames);
    if check_reparam
        Algorithms=["w","a","rep","repA","repAPen","repPen","man","manPen","EM","wMan"];
    else 
        Algorithms=["w","a","man","manPen","EM","wMan"];
    end 
    Algorithms=Algorithms.reshape([],1);
    out_table=[table(Algorithms) out_table];
    
    fig = uifigure("Position",[100,100,760,360]);
    out_table=uitable(fig,"Data",out_table,"Position",[20,20,720,320]);
end 
function res = adapt_res(res,tolgradnorm,minstepsize,maxiter)
    %PROBLEM! Minstepsize is allowed to be under minstepsize 
    %--> Dont know if algorithm would have stopped there
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
        disp(iters)
        disp(k)
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
        disp(iters)
        res.singular=false;
    end 
    k=min([maxiter+1,iters]);
    s=size(names,1);
    names=names(1:s-1);
    for i=1:numel(names)
        name_i=string(names(i));
        if ~isempty(res.(name_i))
            res.(name_i)=res.(name_i)(1:k);
        end 
    end 
end