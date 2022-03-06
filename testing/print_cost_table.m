function print_cost_table
    %General script to print latex cost tables

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
    
    n_results=n_sets*3;
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
    
    iters_man=zeros(n_results,1);
    iters_EM=zeros(n_results,1);
    grad_evals_man=zeros(n_results,1);
    
    relpath_store_basic = "./Results/results";
    count_total=0;
    count_same=0;
    for i=1:size(Ns,2)
        Nstr = "N"+int2str(Ns(i));
        for j=1:size(Ks,2)
            Kstr="K"+int2str(Ks(j));
            %disp(Kstr +"und"+ Dstr)
            for k=1:size(Ds,2)
                Dstr="D"+int2str(Ds(k));
                count=0;
                for l=1:size(cs,2)
                    if cs(l)==0.2
                        cstr="c"+"02";
                    else
                        cstr="c"+int2str(cs(l));
                    end
                    for m=1:n_sets
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
                        
                        if abs(res_w.costs(end)-res_wMan.costs(end))<1e-4
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
                        iters_EM(count)=numel(res_EM.costs);
                        grad_evals_man(count)=sum(res_man.costevals);
                    end
                    

                    %disp("c="+cstr+","+"N="+Nstr+"&"+mean(costs_man)+"\pm"+var(costs_man)+"&"+...
                        %mean(costs_EM)+"\pm"+var(costs_EM)+"&"+mean(times_man)+"\pm"+...
                        %var(times_man)+"&"+mean(times_EM)+"\pm"+var(times_EM)+"\\")
                    %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(times_man),2)+"$\pm$"+...
                     %  round(var(times_man),2)+"&"+round(mean(times_EM),2)+"$\pm$"+round(var(times_EM),2)+"\\")
                    %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(iters_EM),2)+"$\pm$"+...
                     %  round(var(iters_EM),2)+"&"+round(mean(iters_man),2)+"$\pm$"+round(var(iters_man),2)+...
                      % "&"+round(mean(grad_evals_man),2)+"$\pm$"+round(var(grad_evals_man),2)+"\\")
                   
                    %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(costs_man),2)+"$\pm$"+round(var(costs_man),2)+"&"+...
                     %   round(mean(costs_EM),2)+"$\pm$"+round(var(costs_EM),2)+"\\")
                     %disp("D="+Dstr+","+"N="+Nstr+"&"+round(mean(costs_EM),2)+"$\pm$"+round(var(costs_EM),2)+"&"+...
                      %  round(mean(costs_man),2)+"$\pm$"+round(var(costs_man),2)+"\\")
                     
                    %disp("D="+Dstr+","+"N="+Nstr+"&"+mean(iters_EM)+"&"+mean(iters_man)+"&"+...
                     %    mean(grad_evals_man)+"\\");
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
                iters_EM=iters_EM(1:count);
                grad_evals_man=grad_evals_man(1:count);

                %Only take costs not NaN
                times_EM=times_EM(~isnan(times_EM));
                costs_man=costs_man(~isnan(costs_man));
                costs_manPen=costs_manPen(~isnan(costs_manPen));
                costs_wMan=costs_wMan(~isnan(costs_wMan));
                disp("N="+int2str(Ns(i))+","+"K="+int2str(Ks(j))+","+"D="+int2str(Ds(k))+"&"+mean(costs_man)+"$\pm$"+var(costs_man)+"&"+...
                    mean(costs_EM)+"$\pm$"+var(costs_EM)+"&"+mean(times_man)+"$\pm$"+...
                    var(times_man)+"&"+mean(times_EM)+"$\pm$"+var(times_EM)+"\\")
            end 
        end
    end 
    disp(count_same/count_total)
    
   
end 
function res = adapt_res(res,tolgradnorm,minstepsize,maxiter)
    %Function problematic unless same stopping criteria for minstepsize -->
    %Did not actually use it
    
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