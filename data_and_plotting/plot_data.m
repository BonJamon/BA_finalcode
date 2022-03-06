function [] = plot_data(data,mu,Sigma,w)
    if class(Sigma)=="cell"
        D=size(Sigma{1},1);
        K=size(Sigma,1);
        Sigma1 = zeros(D,D,K);
        for i=1:K
            Sigma1(:,:,i)=Sigma{i};
        end 
        Sigma=Sigma1;
    else 
        K=size(Sigma,3);
    end
    if size(w,2)==(K-1)
        w = [w,0];
        p = exp(w-max(w)) ./ sum(exp(w-max(w)));
    else
        p=w;
    end     
 
    figure;
    scatter(data(:,1),data(:,2));
    hold on
    x1min = min(data(:,1));
    x1max = max(data(:,1));
    x2min = min(data(:,2));
    x2max = max(data(:,2));
    [X1,X2] = meshgrid(x1min:0.1:x1max,x2min:0.1:x2max);
    len1 = size(X1,1);
    len2 = size(X1,2);
    pos =[reshape(X1,len1*len2,1),reshape(X2,len1*len2,1)];
    colours = ["r","g","b","c","k","m"];
    for i=1:K
        vals = mvnpdf(pos,mu(i,1:2),Sigma(1:2,1:2,i));
        ymax = max(vals);
        vals = reshape(vals,len1,len2);
        contour(X1,X2,vals,[0.2*ymax,0.4*ymax,0.6*ymax,0.8*ymax,ymax],colours(i));
    end 
    