function comp_wMan_Man
    %General script to compare results vor lbfgs_w, lbfgs_wMan, lbfgs_man

    %settings
    tolgradnorm=1e-6;
    minstepsize=1e-12;
    maxiter=1500;
    
    %init
    Ns=[250,2500,25000];
    Ks=[2,5];
    Ds=[2,5,20];
    cs=[0.2,1,5];
    n_sets=5;
    
    n_results=n_sets*size(Ks,2)*size(cs,2);
    n_resultsAll=n_sets*size(Ks,2)*size(cs,2)*size(Ds,2)*size(Ns,2);
    costs_w=zeros(n_results,1);
    times_w=zeros(n_results,1);
    costs_a=zeros(n_results,1);
    times_a=zeros(n_results,1);

    costs_man=zeros(n_results,1);
    times_man=zeros(n_results,1);
    costs_manPen=zeros(n_results,1);
    times_manPen=zeros(n_results,1);
    costs_EM=zeros(n_results,1);
    times_EM=zeros(n_results,1);
    costs_wMan=zeros(n_results,1);
    times_wMan=zeros(n_results,1);
    
    costs_wManAll=zeros(n_resultsAll,1);
    costs_manAll=zeros(n_resultsAll,1);
    times_wManAll=zeros(n_resultsAll,1);
    times_manAll=zeros(n_resultsAll,1);
    
    times_repAll=zeros(n_resultsAll,1);
    iters_repAll=zeros(n_resultsAll,1);
    
    iters_man=zeros(n_results,1);
    iters_wMan=zeros(n_results,1);
    iters_w=zeros(n_results,1);
    gradevals_wMan=zeros(n_results,1);
    gradevals_w=zeros(n_results,1);
    gradevals_man=zeros(n_results,1);
    
    relpath_store_basic = "./Results/results";
    count_total=0;
    count_same=0;
    countAll=0;
    for k=1:size(Ds,2)
        Dstr="D"+int2str(Ds(k));
        for i=1:size(Ns,2)
            Nstr = "N"+int2str(Ns(i));
            count=0;
            for l=1:size(cs,2)
                if cs(l)==0.2
                    cstr="c"+"02";
                else
                    cstr="c"+int2str(cs(l));
                end
                for j=1:size(Ks,2)
                    Kstr="K"+int2str(Ks(j));
                    for m=1:n_sets
                        countAll=countAll+1;
                        count_total=count_total+1;
                        relpath_store=relpath_store_basic+"/"+Nstr+"/"+Kstr+"/"+Dstr+"/"+cstr;
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
                            res_EM=load(filename_EM).res_EM;
                            res_EM=adapt_resEM(res_EM,maxiter);
                        end 
                        
                        filename_wMan = relpath_store+"/wMan_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_wMan)
                            res_wMan = load(filename_wMan).res_w_man;
                            res_wMan=adapt_res(res_wMan,tolgradnorm,minstepsize,maxiter);
                        end 
                        filename_rep = relpath_store+"/rep_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_rep)
                            res_rep = load(filename_rep).res_rep;
                            res_rep=adapt_res(res_rep,tolgradnorm,minstepsize,maxiter);
                        end 
                        
                        if abs(res_w.costs(end)-res_wMan.costs(end))<1e-4 && abs(res_wMan.costs(end)-res_man.costs(end))<1e-4&&abs(res_man.costs(end)-res_w.costs(end))<1e-4
                            count_same=count_same+1;
                        end 
                        costs_w(count)=res_w.costs(end);
                        times_w(count)=res_w.times(end);
                        costs_a(count)=res_a.costs(end);
                        times_a(count)=res_a.times(end);
                        costs_man(count)=res_man.costs(end);
                        times_man(count)=res_man.times(end);
                        costs_manPen(count)=res_man_pen.costs(end);
                        times_manPen(count)=res_man_pen.times(end);
                        costs_EM(count)=res_EM.costs(end);
                        %been a mistake in EM code for singular case
                        if isempty(res_EM.times)
                            times_EM(count)=NaN;
                        else 
                            times_EM(count)=res_EM.times(end);
                        end 
                        costs_wMan(count)=res_wMan.costs(end);
                        times_wMan(count)=res_wMan.times(end);
                        
                        iters_man(count)=numel(res_man.costs);
                        iters_wMan(count)=numel(res_wMan.times);
                        iters_w(count)=numel(res_w.costs);
                        gradevals_man(count)=sum(res_man.costevals);
                        gradevals_wMan(count)=sum(res_wMan.costevals);
                        gradevals_w(count)=sum(res_w.costevals);
                        

                        costs_wManAll(countAll)=res_wMan.costs(end)/res_w.costs(end);
                        costs_manAll(countAll)=res_man.costs(end)/res_w.costs(end);
                        times_wManAll(countAll)=res_wMan.times(end);%res_wMan.times(end)/res_w.times(end);
                        times_manAll(countAll)=res_man.times(end);%/res_w.times(end);
                        
                        times_repAll(countAll)=res_rep.times(end)/res_w.times(end);
                        iters_repAll(countAll)=numel(res_rep.costs)/numel(res_w.costs);
                    end
                end
            end 
            %For the case not all results are there
            costs_w=costs_w(1:count);
            times_w=times_w(1:count);
            costs_a=costs_a(1:count);
            times_a=times_a(1:count);   
            costs_man=costs_man(1:count);
            times_man=times_man(1:count);
            costs_manPen=costs_manPen(1:count);
            times_manPen=times_manPen(1:count);
            costs_EM=costs_EM(1:count);
            times_EM=times_EM(1:count);
            costs_wMan=costs_wMan(1:count);
            times_wMan=times_wMan(1:count);
            iters_man=iters_man(1:count);
            iters_wMan=iters_wMan(1:count);
            iters_w=iters_w(1:count);
            gradevals_man=gradevals_man(1:count);
            gradevals_wMan=gradevals_wMan(1:count);
            gradevals_w=gradevals_w(1:count);

            %Only take costs not NaN
            times_EM=times_EM(~isnan(times_EM));
            costs_man=costs_man(~isnan(costs_man));
            costs_manPen=costs_manPen(~isnan(costs_manPen));
            costs_wMan=costs_wMan(~isnan(costs_wMan));

            
            %disp("c="+cstr+","+"N="+Nstr+"&"+mean(costs_man)+"\pm"+var(costs_man)+"&"+...
                %mean(costs_EM)+"\pm"+var(costs_EM)+"&"+mean(times_man)+"\pm"+...
                %var(times_man)+"&"+mean(times_EM)+"\pm"+var(times_EM)+"\\")
            %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(times_wMan),2)+"$\pm$"+...
             %  round(var(times_wMan),2)+"&"+round(mean(times_man),2)+"$\pm$"+round(var(times_man),2)+"\\")
            %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(gradevals_wMan),2)+"\\")
             %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(iters_wMan),2)+"$\pm$"+...
             %  round(var(iters_wMan),2)+"&"+round(mean(iters_man),2)+"$\pm$"+round(var(iters_man),2)+...
              % "&"+round(mean(gradevals_wMan),2)+"$\pm$"+round(var(gradevals_wMan),2)+...
               %"&"+round(mean(gradevals_man),2)+"$\pm$"+round(var(gradevals_man),2)+"\\")

             %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(costs_wMan),2)+"$\pm$"+round(var(costs_wMan),2)+"&"+...
              %  round(mean(costs_man),2)+"$\pm$"+round(var(costs_man),2)+"\\")
              
             %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(costs_w),2)+"$\pm$"+round(var(costs_w),2)+...
              %   "&"+round(mean(costs_wMan),2)+"$\pm$"+round(var(costs_wMan),2)+"&"+...
               %  round(mean(costs_man),2)+"$\pm$"+round(var(costs_man),2)+"\\")
            
            %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(gradevals_w),2)+"$\pm$"+round(var(gradevals_w),2)+...
             %    "&"+round(mean(gradevals_wMan),2)+"$\pm$"+round(var(gradevals_wMan),2)+"&"+...
              %  round(mean(gradevals_man),2)+"$\pm$"+round(var(gradevals_man),2)+"\\")
              

            %disp("D="+Dstr+","+"N="+Nstr+"&"+mean(iters_EM)+"&"+mean(iters_man)+"&"+...
             %    mean(grad_evals_man)+"\\");
        end
        disp("D="+Dstr+"&"+round(mean(gradevals_w),2)+"$\pm$"+...
             round(var(gradevals_w),2)+"&"+round(mean(gradevals_wMan),2)+"$\pm$"+round(var(gradevals_wMan),2)+"\\")
        
    end 
    costs_wManAll=costs_wManAll(1:countAll);
    costs_manAll=costs_manAll(1:countAll);
    times_wManAll=times_wManAll(1:countAll);
    times_manAll=times_manAll(1:countAll);
    costs_wManAll=costs_wManAll(~isnan(costs_wManAll));
    costs_manAll=costs_manAll(~isnan(costs_manAll));
    
    %format short
    %disp(mean(costs_wManAll))
    %disp(mean(times_wManAll))
    %disp(mean(costs_manAll))
    %disp(mean(times_manAll))
    %disp(count_same/count_total)
    %disp(mean(times_repAll))
    %disp(mean(iters_repAll))
    disp(sum(times_wManAll)/sum(times_manAll))
    
   
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