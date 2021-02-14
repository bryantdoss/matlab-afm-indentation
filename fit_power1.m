function [fitobj, gof] = fit_power1(x, y)
    if length(x) ~= length(y)
        disp('Fit problem: lengths not the same')
        return
    end
    
    % Fits y = A*x^B
    
    % Linear least squares
%     n = length(x);
%     fitobj.b = (n * sum(log(x).*log(y)) - sum(log(x))*sum(log(y))) / (n*sum(log(x).^2) - sum(log(x))^2);
%     fitobj.a = exp((sum(log(y)) - fitobj.b*sum(log(x))) / n);
%     gof = nan;    
    
%     % Non-linear numeric least squares
%     
%     % Scale x and y
%     scale = 1 / max(y);
%     xs = x * scale;
%     ys = y * scale;
%     
%     % Non-linear fitting function
%     %f = @(xf, xd) xf(2).*xd.^xf(1);
%     f = @(xf, xd) xf(2).*xd.^xf(1) + xf(3);
%     
%     % Initial power law value from linear least squares
%     n = length(x);
%     b = (n * sum(log(x).*log(y)) - sum(log(x))*sum(log(y))) / (n*sum(log(x).^2) - sum(log(x))^2);
%     b = 1.5;
%     
%     % Perform the fit
%     %fits = lsqcurvefit(f, [b 0], xs, ys, [0 0], [5 Inf], optimset('TolFun', 1e-9, 'Display', 'off', 'Algorithm','trust-region-reflective'));
%     fits = lsqcurvefit(f, [b 0 0], xs, ys, [0 0 -Inf], [5 Inf Inf], optimset('TolFun', 1e-9, 'Display', 'off', 'Algorithm','trust-region-reflective'));
%     
%     % Rescale the fits
%     fitobj.b = fits(1);
%     fitobj.a = fits(2) * (1/scale) / (1/scale)^fits(1);
%     %fitobj.c = fits(3) * (1/scale);
%     %gof = nan;
    
    if x(1) < 10*eps % Improve robustness of fit
        x(1) = [];
        y(1) = [];
    end

    % Use fit_poly1
    yf = log(y);
    xf = log(x);
    [fit1,gof1] = fit_poly1(xf,yf);
    fitobj.b = fit1.p1;
    fitobj.a = exp(fit1.p2);
    gof.b_se = gof1.p1_se;
    gof.a_se = gof1.p2_se*exp(fit1.p2);
    
%     figure
%     plot(x,y,'linewidth',2)
%     hold on
%     plot(x, fitobj.a.*x.^fitobj.b,'k');
%     figure
%     plot(x, y - (fitobj.a.*x.^fitobj.b));
end
