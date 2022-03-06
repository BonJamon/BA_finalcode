function [out_table,plots] = evaluate_results(N,K,D,c)
    %Function evaluating results for data simulated according to (N,K,D,c)
    %Produces: table containing cost, time, iterations, costevaluations
    %taken AND Plots for Median of results for different algorithms showing
    %origin of improvement

    check_reparam=true;
    Nstr="N"+int2str(N);
    Kstr="K"+int2str(K);
    Dstr="D"+int2str(D);
    if c==0.2
        cstr="c"+"02";
    else 
        cstr="c"+int2str(c);
    end
    path="./Results/results/"+Nstr+"/"+Kstr+"/"+Dstr+"/"+cstr;
    partstr=Nstr+Kstr+Dstr+cstr;
    
    if check_reparam
        algorithms=["w","a","rep","repA","repPen","repAPen","man","manPen","EM","wMan"];
    else 
        algorithms=["w","a","man","manPen","EM","wMan"];
    end
    mean_indices=zeros(1,size(algorithms,2));
    
    %Get the table
    out_table=zeros(size(algorithms,2),8);
    for i=1:size(algorithms,2)
        costs=zeros(1,5);
        times=zeros(1,5);
        iters=zeros(1,5);
        costevals=zeros(1,5);
        for j=1:5
            filename=path+"/"+algorithms(i)+"_"+partstr+"_"+int2str(j);
            res=load(filename);
            res=struct2cell(res);
            res=res{1};
            costs(j)=res.costs(end);
            times(j)=res.times(end);
            iters(j)=numel(res.costs);
            if i~=size(algorithms,2)-1
                costevals(j)=sum(res.costevals);
            end
        end 

       
        out_table(i,:)=[mean(costs),std(costs),mean(times),std(times),...
            mean(costevals),std(costevals),mean(iters),std(iters)];
        %Problem: When all indices are NaN-->rn just take mean_index=1
        mean_index=get_mean_index(costs);
        mean_indices(i)=mean_index;
    end 
    varnames=["Mean Cost","std ofCost","Mean Time",...
        "std of Time","Mean Costevaluations","std of Costevaluations",...
"Mean Number of Iterations","std of Number of Iterations"];
    out_table=array2table(out_table,"VariableNames",varnames);
    Algorithm=reshape(algorithms,size(algorithms,2),1);
    out_table=[table(Algorithm) out_table];
    
    fig = uifigure("Position",[100,100,760,760]);
    out_table=uitable(fig,"Data",out_table,"Position",[20,20,720,720]);
    
    %Get the plots
    results=cell(1,size(algorithms,2));
    for i=1:size(algorithms,2)
        filename=path+"/"+algorithms(i)+"_"+partstr+"_"+int2str(mean_indices(i));
        res=load(filename);
        res=struct2cell(res);
        res=res{1};
        results{i}=res;
    end 
    
    costs_w=results{1}.costs;
    times_w=results{1}.times;
    costs_a=results{2}.costs;
    times_a=results{2}.times;
    if check_reparam
        costs_rep=results{3}.costs;
        costs_repA=results{4}.costs;
        costs_repPen=results{5}.costs;
        costs_repAPen=results{6}.costs;
    end 
    costs_man=results{7}.costs;
    costs_manPen=results{8}.costs;
    costs_EM=results{9}.costs;
    costs_wMan=results{10}.costs;
    
    iter_w = size(costs_w,2);
    iterations_w = linspace(0,iter_w-1,iter_w);
    iter_a = size(costs_a,2);
    iterations_a = linspace(0,iter_a-1,iter_a);
    if check_reparam
        iter_rep = size(costs_rep,2);
        iterations_rep = linspace(0,iter_rep-1,iter_rep);
        iter_repA = size(costs_repA,2);
        iterations_repA = linspace(0,iter_repA-1, iter_repA);
        iter_repPen = size(costs_repPen,2);
        iterations_repPen=linspace(0,iter_repPen-1,iter_repPen);
        iter_repAPen = size(costs_repAPen,2);
        iterations_repAPen=linspace(0,iter_repAPen-1,iter_repAPen);
    end 
    iter_man = size(costs_man,2);
    iterations_man = linspace(0,iter_man-1,iter_man);
    iter_manPen = size(costs_manPen,2);
    iterations_manPen = linspace(0,iter_manPen-1,iter_manPen);
    iter_EM=size(costs_EM,1);
    iterations_EM = linspace(0,iter_EM-1,iter_EM);
    iter_wMan=size(costs_wMan,2);
    iterations_wMan=linspace(0,iter_wMan-1,iter_wMan);
    costs_EM=reshape(costs_EM,1,iter_EM);
    
    if check_reparam
        plots(1)=figure;
        %plot(iterations_w,costs_w,"c",iterations_a,costs_a,"m",...
         %   iterations_rep, costs_rep,"b",iterations_repA,...
          %  costs_repA,"r",iterations_EM,costs_EM,"k");
         hold on
         plot(times_w,costs_w,"c",times_a, costs_a,"b");
         plot(times_w(end),costs_w(end),"c*");
         plot(times_a(end),costs_a(end),"b*");
         
        xlabel("Iteration");
        ylabel("Cost");
        title("Cost vs Time: Influence of multinomial manifold");
        %legend("LBFGS using w","LBFGS using a","Reparametrisized LBFGS using w",...
           % "Reparametrisized LBFGS using a","EM");
        legend("LBFGS using euclidean mixing weights", "LBFGS using multinomial manifold");
        hold off
        plots(2)=figure;
        plot(iterations_w,costs_w,"c",iterations_rep,costs_rep,"b",iterations_repPen,costs_repPen,"g",...
            iterations_a,costs_a,"m",iterations_repA,costs_repA,"r",...
            iterations_repAPen,costs_repAPen,"y",iterations_EM,costs_EM,"k");
        xlabel("Iteration");
        ylabel("Cost");
        title("Cost vs Iteration: Influence of Penalizer");
        legend("LBFGS using w","Reparametrisized LBFGS using w",...
            "Reparametrisized LBFGS penalized","LBFGS using a","Reparametrizised LBFGS using a",...
            "Reparametrisized LBFGS using a penalized","EM");

        plots(3)=figure;
        plot(iterations_w,costs_w,"c",iterations_rep,costs_rep,"b",iterations_repPen,costs_repPen,"g",...
        iterations_man, costs_man,"m",iterations_manPen, costs_manPen,"r",iterations_EM,costs_EM,"k",...
        iterations_wMan,costs_wMan);
        xlabel("Iteration");
        ylabel("Cost");
        title("Cost vs Iteration: Influence of Manifold");
        legend("LBFGS using w","Reparametrisized LBFGS using w",...
            "Reparametrisized LBFGS penalized","Manifold LBFGS","Manifold LBFGS penalized","EM",...
            "Manifold LBFGS on original formulation");
    else
        plots(1)=figure;
        plot(iterations_w,costs_w,"c",iterations_man, costs_man,"g",iterations_EM,costs_EM,"k",...
        iterations_wMan,costs_wMan,"r");
        xlabel("Iteration");
        ylabel("Cost");
        title("Cost vs Iteration: Influence of Manifold");
        legend("LBFGS using w","Manifold LBFGS","EM",...
            "Manifold LBFGS on original formulation");
        
        %TODO: Input costs aMan when implemented
        plots(2)=figure;
        plot(iterations_w,costs_w,"c",iterations_a,costs_a,"b",...
            iterations_wMan, costs_wMan,"r",iterations_man,...
            costs_man,"g",iterations_EM,costs_EM,"k");
        xlabel("Iteration");
        ylabel("Cost");
        title("Cost vs Iteration: Influence of Reparametrisation");
        legend("LBFGS using w","LBFGS using a","orginal Manifold LBFGS",...
            "reparametrisized Manifold LBFGS","EM");

        plots(3)=figure;
        plot(iterations_w,costs_w,"c",iterations_wMan,costs_wMan,"r",iterations_man,costs_man,"g",...
            iterations_manPen,costs_manPen,"y",iterations_EM,costs_EM,"k");
        xlabel("Iteration");
        ylabel("Cost");
        title("Cost vs Iteration: Influence of Penalizer");
        legend("LBFGS using w","original Manifold LBFGS",...
            "Reparametrisized Manifold LBFGS","penalized reparametrisized Manifold LBFGS","EM");
 
    end
end

function mean_index = get_mean_index(costs)
    %Gives the mean index for all cost values not NaN
    costsnotnan=costs(~isnan(costs));
    disp(costs)
    %Problem: When costsnotnan%2=0 --> Takes mean of middle two values
    mediancost=median(costsnotnan(1,:));
    bool1=abs(costs(1,:)-mediancost)==min(abs(costs(1,:)-mediancost));
    bool2=~isnan(costs(1,:));
    mean_index=find(bool1&bool2,1);
    if isempty(mean_index)
        mean_index=1;
    end 
end 
    
    