function [fitobj, gof] = fit_poly1(x, y)
    if length(x) ~= length(y)
        disp('Fit problem: lengths not the same')
        return
    end
    n = length(x);
    fitobj.p1 = (mean(x.*y) - mean(x)*mean(y)) / (mean(x.^2) - mean(x)^2);
    fitobj.p2 = mean(y) - fitobj.p1*mean(x);
    
    fit = x * fitobj.p1 + fitobj.p2;
    
%     if x(1) > 3.5e-6
%     plot(x, fit, x, y)
%     return
%     end
    
    ssres = sum((y - fit).^2);
    sstot = sum((y-mean(y)).^2);
    fitobj.r2 = 1 - ssres/sstot;

    % Calculate 99% confidence interval error for p1 (the slope)
    % s = sqrt((1/(n-2)*sum((fit-y).^2))/(sum((x-mean(x)).^2)));
    % gof.p1 = s*tquant(n-2,1-.01/2);
    
    % Calculate standard error
%     gof.sse = sum((fit-y).^2);
%     gof.mse = gof.sse / (length(x)-2);
%     gof.p1_se = sqrt(gof.mse / sum((x - mean(x)).^2));
%     gof.p2_se = sqrt(gof.mse) * sqrt(1/length(x) - mean(x)^2 / sum((x - mean(x)).^2));
%     
%     gof.p1_sd = gof.p1_se * sqrt(n);
%     gof.p2_sd = gof.p2_se * sqrt(n);
%     
%     gof.pcc = sum((x-mean(x)).*(y-mean(y)))/sqrt(sum((x-mean(x)).^2)*sum((y-mean(y)).^2));
%     
%     gof.pcc = sqrt(sum((fit-mean(y)).^2)/sum((y-mean(y)).^2));

    gof.ssxx = sum((x-mean(x)).^2);
    gof.ssyy = sum((y-mean(y)).^2);
    gof.ssxy = sum((x-mean(x)).*(y-mean(y)));
    gof.s = sqrt((gof.ssyy-fitobj.p1*gof.ssxy)/(n-2));
    
    gof.p2_se = gof.s * sqrt(1/n+mean(x)^2/gof.ssxx);
    gof.p1_se = gof.s / sqrt(gof.ssxx);

end
