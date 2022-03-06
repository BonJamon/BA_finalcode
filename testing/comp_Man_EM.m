function fig = comp_man_EM()
    %General script to compare man and EM results
    
    %settings
    tolgradnorm=1e-6;
    minstepsize=1e-12;
    maxiter=1500;
    relational=false;
    comp_iters=false;
    
    %init
    count=0;
    Ns=[250,2500,25000];
    Ks=[2,5];
    Ds=[2,5,20];
    cs=[0.2,1,5];
    n_sets=5;
    
    n_results=size(Ns,2)*size(Ks,2)*size(Ds,2)*size(cs,2)*n_sets;

    costs_man=zeros(n_results,1);
    times_man=zeros(n_results,1);
    iters_man=zeros(n_results,1);
    costs_EM=zeros(n_results,1);
    times_EM=zeros(n_results,1);
    iters_EM=zeros(n_results,1);
    cost_evals_man=[];

    
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
                        
                        filename_EM = relpath_store+"/EM_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_EM)
                            count=count+1;
                            disp(filename_EM)
                            res_EM=load(filename_EM).res_EM;
                            res_EM=adapt_resEM(res_EM,maxiter);
                        else 
                            %Been a mistake, where EM couldnt get saved bc
                            %file too big
                            continue;
                        end 
                        
                        
                        filename_man = relpath_store+"/man_"+Nstr+Kstr+Dstr+cstr+"_"+int2str(m)+".mat";
                        if isfile(filename_man)
                            res_man=load(filename_man).res_man;
                            res_man=adapt_res(res_man,tolgradnorm,minstepsize,maxiter);
                        end 
                        
                        if relational
                            costs_man(count)=1;
                            times_man(count)=1;
                            iters_man(count)=1;
                            costs_EM(count)=res_EM.costs(end)/res_man.costs(end);
                            iters_EM(count)=numel(res_EM.costs)/numel(res_man.costs);
                            %been a mistake in EM code for singular case
                            if isempty(res_EM.times)
                                times_EM(count)=NaN;
                            else 
                                times_EM(count)=res_EM.times(end)/res_man.times(end);
                            end
                            cost_evals_man(count)=sum(res_man.costevals);
                        else 
                            costs_man(count)=res_man.costs(end);
                            times_man(count)=res_man.times(end);
                            iters_man(count)=numel(res_man.costs);
                            iters_EM(count)=numel(res_EM.costs);
                            costs_EM(count)=res_EM.costs(end);
                            if isempty(res_EM.times)
                                times_EM(count)=NaN;
                            else 
                                times_EM(count)=res_EM.times(end);
                            end
                            cost_evals_man(count)=sum(res_man.costevals);
                            %[cost_evals_man, res_man.costevals];
                        end
                    end
                end
            end 
        end
    end 
    costs_man=costs_man(1:count);
    times_man=times_man(1:count);
    iters_man=iters_man(1:count);
    costs_EM=costs_EM(1:count);
    times_EM=times_EM(1:count);
    iters_EM=iters_EM(1:count);
    cost_evals_man=cost_evals_man(1:count);
    
    times_EM=times_EM(~isnan(times_EM));
    costs_EM=costs_EM(~isnan(costs_EM));
    costs_man=costs_man(~isnan(costs_man));
    disp("costevals")
    disp(mean(cost_evals_man))
    disp(var(cost_evals_man))
    disp("Iters man")
    disp(mean(iters_man))
    disp(var(iters_man))
    disp("times man")
    disp(mean(times_man))
    disp(var(times_man))
    disp("iters EM")
    disp(mean(iters_EM))
    disp(var(iters_EM))
    disp("time EM")
    disp(mean(times_EM))
    disp(var(times_EM))
    

    
    
    if ~comp_iters
        stats(1,:)=[mean(costs_man),std(costs_man),mean(times_man),std(times_man)];
        stats(2,:)=[mean(costs_EM),std(costs_EM),mean(times_EM),std(times_EM)];
    else
        stats(1,:)=[mean(costs_man),std(costs_man),mean(iters_man),std(iters_man)];
        stats(2,:)=[mean(costs_EM),std(costs_EM),mean(iters_EM),std(iters_EM)];
    end
    if ~comp_iters
        if relational
            varNames=["Ratio to man cost","std of Ratio to man cost","Mean Ratio to man time","std of Ratio to man time"];
        else
            varNames=["Mean Cost","Std Cost","Mean Time","Std Time"];
        end
    else 
        if relational
            varNames=["Ratio to man cost","std of Ratio to man cost","Mean ratio to man iters","std of Ratio to man iters"];
        else 
            varNames=["Mean Cost","Std Cost","Mean Iters","Std Iters"];
        end 
    end 
    out_table=array2table(stats,"VariableNames",varNames);
    Algorithms=["man","EM"];
    Algorithms=Algorithms.reshape([],1);
    out_table=[table(Algorithms) out_table];
    
    fig = uifigure("Position",[100,100,760,360]);
    out_table=uitable(fig,"Data",out_table,"Position",[20,20,720,320]);
end 
function res = adapt_res(res,tolgradnorm,minstepsize,maxiter)
    %Find index manhere it manould have stopped
    t=find(res.gradnorms<tolgradnorm);
    if isempty(t)
        t=maxiter;
    else 
        t=t(1);
    end 
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
        res.singular=flase;
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