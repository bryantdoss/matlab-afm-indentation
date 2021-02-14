% dd_fit.m
%
% Depth-dependent force curve fitting
% Input parameters: 
%  curve: the force curve
%  contact: the contact point
%  r: the tip radius
%  theta: axisymmetric semi-vertical tip angle  
%  nu: Poisson ratio
%  binning: bins for depth-dependent analysis
%  type: the indenter probe geometry type, supported types:
%       parabolic: Hertz model
%       conical: Sneddon model
%       hyperbolic: a hyperbola
%       briscoe: Briscoe's blunted cone model
%       spherocone: sphrical tip, transitions to cone continuosly
%       cylindrical: a cylinder
%       sphere: a sphere
%  min_ind: Minimum indentation depth
%  max_ind: Maxiumum indentation depth
% Output parameters
%  E: the depth-dependent Young's modulus
%  X: the indentation value so plot(X,E) works
%  r2: the r^2 value from the fit
% Update 2020: now the zero force is chosen from an "average" of baseline,
% also power-law fitting on contact model excludes very small 
% force/indentation values to make it more robust

function [E, X, r2, offset] = dd_fit(curve, contact, r, theta, nu, binning, type, min_ind, max_ind, doplot, h)

if nargin < 9
    min_ind = 0;
    max_ind = 42;
    doplot = 'no';
    h = 0;
end

if nargin < 10
    doplot = 'no';
    h = 0;
end

% if length(curve) > 1 && binning == 1
%     X = [];
%     r2 = [];
%     for i=1:length(curve)
%         E(i) = dd_fit(curve(i), contact(i), r, theta, nu, binning, type, min_ind, max_ind, doplot, h);
%     end
%     return
% end

if size(curve.ext,1) < 20 || isnan(contact)
    E = NaN;
    X = nan;
    r2 = nan;
    offset = nan;
    return
end

% Pre-process the force curve, truncate the baseline and shift
% everything to zero
index = get_index(curve.ext(:,1),contact);
index2 = get_index(curve.ext(:,1),contact + max_ind);
ext_tmp = curve.ext(1:index,:);
zeroF = mean(ext_tmp(:,2));
ext = curve.ext(index:index2, :);
ext(:,1) = ext(:,1) - min(ext(:,1)) + eps;
% ext(:,2) = ext(:,2) - min(ext(:,2)) + eps;
% ext(:,2) = ext(:,2) - ext(1,2) + eps;
ext(:,2) = ext(:,2) - zeroF;
for i=size(ext,1)-1:-1:1
    if ext(i,2) < 0
        min_ind = ext(i+1,1)+eps;
        break
    end
    if ext(i,1) < min_ind
        break
    end
end
fc = curve;
fc.ext = ext;
% if (size(fc.ext,1) < 2)
%     fc.ext = curve.ext(end-4:end,:);
% end
% index3 = get_index(fc.ext(:,1), min_ind);
% ext2(:,1) = ext(index3:end,1);
% ext2(:,2) = ext(index3:end,2);
% fc.ext = ext2;
% ext = ext2;

% Initialize storage for contact radii and lambda terms
as = [];
lams = [];

if length(ext) < 3
    E = NaN;
    X = nan;
    r2 = NaN;
    offset = nan;
    return
end

% Iterate over indentation depths
for i=1:length(ext)
    d = ext(i,1);
    if strcmp(type,'parabolic')
        a = sqrt(r*d);
        lam = (4/3) * sqrt(r) * d^(3/2);
        if h ~= 0
            % Garcia & Garcia correction
            lam = lam .* (1 + 1.133*sqrt(r*d)/h + 1.497*r*d/h^2 + 1.469*(r*d)^(3/2)/h^3 + 0.755*r^2*d^2/h^4);
        end        
    elseif strcmp(type,'conical')
        a = 2*d/(pi*cot(theta));
        lam = (1/2) * a^2 * pi * cot(theta);
    elseif strcmp(type,'hyperbolic')
        f = @(a) real(d - (a*cot(theta)/2) * (pi/2 + atan(a/(2*r*cot(theta)) - r*cot(theta)/(2*a))));
        a = sqrt(r*d);
        a = fzero(f,a);
        z = r*cot(theta)/a;
        lam = (a^3/r)*(z^2 + (z/2)*(1-z^2)*(pi/2 + atan(1/(2*z)-(z/2))));
    elseif strcmp(type,'briscoe')
        b = r * cos(theta);
        m = 1/2;
        n = 1;
        f = @(a) real((ext(i,1) + (a/r)*((a^2-b^2)^(1/2)-a) - ...
                 (n*a/tan(theta))*((pi/2)-asin(b/a))));
        a = sqrt(r*ext(i,1));
        a = fzero(f,a);
        lam = 2 * (a*d - (m*a^2/tan(theta))*(pi/2 - asin(b/a)) ...
            - a^3/(3*r) + ((a^2-b^2)^(1/2))*(m*b/tan(theta)+(a^2-b^2)/(3*r)));
    elseif strcmp(type,'rico')
        b = r * cos(theta);
        m = sqrt(2)/pi;
        n = 2^(3/2)/pi;
        f = @(a) real((ext(i,1) + (a/r)*((a^2-b^2)^(1/2)-a) - ...
                 (n*a/tan(theta))*((pi/2)-asin(b/a))));
        a = sqrt(r*ext(i,1));
        a = fzero(f,a);
        lam = 2 * (a*d - (m*a^2/tan(theta))*(pi/2 - asin(b/a)) ...
            - a^3/(3*r) + ((a^2-b^2)^(1/2))*(m*b/tan(theta)+(a^2-b^2)/(3*r)));
    elseif strcmp(type,'spherocone')
        b = r*cot(theta) / sqrt(1+cot(theta)^2);
        b = r*cos(theta);
        f = @(a) real(ext(i,1) - ((1 <= b/a) * (1/2)*a*log((r+a)/(r-a)) ...
            + (1 > b/a) * (a*log((r+a)/(sqrt(r^2-b^2)+a*sqrt(1-(b/a)^2))) + a*acos(b/a)*cot(theta))));
        a = sqrt(r*ext(i,1)+1e-9);
        a = fzero(f,a);
        if 1 <= b/a
            lam = ((1/2)*(a^2+r^2)*log((r+a)/(r-a))-r*a);
        elseif 1 > b/a
            lam = a^2*cot(theta)*acos(b/a) + b*cot(theta)*sqrt(a^2-b^2) - a*r + sqrt((r^2-b^2)*(a^2-b^2)) + a^2*log((r+a)/(sqrt(r^2-b^2)+sqrt(a^2-b^2))) - r^2/2*log((a^2*r^2-(b^2-sqrt((r^2-b^2)*(a^2-b^2)))^2)/(b^2*(r+a)^2));
        end
    elseif strcmp(type,'cylindrical')
        a = r;
        lam = 2 * a * d;
    elseif strcmp(type, 'parabocone')
        b = r * cot(theta);
        f = @(a) real(ext(i,1) - ((1 < b/a) * (a^2/r) ...
            + (1 > b/a) * (a^2/r * (1-sqrt(1-b^2/a^2)) + a*acos(b/a)*cot(theta))));
        a = sqrt(r*ext(i,1));
        a = fzero(f,a);
        lam = (1<=b/a) * (2 * a * (d - a^2/(3*r))) ...
            + (1> b/a) * (2 * a * (d - a^2/(3*r)*(1-sqrt(1-b^2/a^2)+b^2/a^2*sqrt(1-b^2/a^2)) + 1/2*a*(b/a*sqrt(1-b^2/a^2)-acos(b/a))*cot(theta)));
    elseif strcmp(type, 'dimitriadis-bonded')
        height = 4e-6;
        a = sqrt(r*d);
        a0 = -(1.2876 - 1.4678*nu + 1.3442*nu^2)/(1-nu);
        b0 =  (0.6387 - 1.0277*nu + 1.5164*nu^2)/(1-nu);
        c = sqrt(r*d)/height;
        lam = (4/3) * sqrt(r) * d^(3/2) * (1 - 2*a0*c/pi + 4*a0^2*c^2/(pi^2) - 8/(pi^3)*(a0^3+4*pi^2*b0/15)*c^3 + 16*a0/(pi^4)*(a0^3+3*pi^2*b0/5)*c^4);
    elseif strcmp(type, 'dimitriadis-nonbonded')
        height = 4e-6;
        a = sqrt(r*d);
        a0 = -0.347 * (3-2*nu)/(1-nu);
        b0 =  0.056 * (5-2*nu)/(1-nu);
        c = sqrt(r*d)/height;
        if c > 1
            disp('warning: chi > 1')
        end
        lam = (4/3) * sqrt(r) * d^(3/2) * (1 - 2*a0*c/pi + 4*a0^2*c^2/(pi^2) - 8/(pi^3)*(a0^3+4*pi^2*b0/15)*c^3 + 16*a0/(pi^4)*(a0^3+3*pi^2*b0/5)*c^4);
    end
    as = [as real(a)];
    lams = [lams real(lam)];
    
end
max_ind = ext(end,1);
nbins = max(floor(max_ind/binning - min_ind/binning),1);

for i=1:nbins
    % Find first point, last point
    first = get_index(fc.ext(:,1), binning*(i-1)+min_ind);
    last = get_index(fc.ext(:,1), binning*i+min_ind);
    % Improve robustness of power-law fit, small numbers hurt it
    firstpower = get_index(fc.ext(:,1), max(fc.ext(first,1), fc.ext(last,1)/6));
    % Fit the model data to a power-law
    [fitobj,gof] = fit_power1(fc.ext(firstpower:last,1), lams(firstpower:last)');
    A = fitobj.a;
    B = fitobj.b;
    % Linearize the experimental data
    ext_lin = ext(first:last,:);
    ext_lin(:,2) = ext_lin(:,2).^(1/B);
    %ext_lin(:,2) = ext_lin(:,2) - min(ext_lin(:,2));
    % Fit linearized curve to a y=mx+b
    [fitobj2,gof2] = fit_poly1(ext_lin(:,1), ext_lin(:,2));
    C = fitobj2.p1;
    D = fitobj2.p2;
    % Calculate the Young's modulus
    E(i) = C^(B) * (1-nu^2) / A;
    offset(i) = -D / C;
    r2(i) = fitobj2.r2;    
    
    if strcmp(doplot,'yes')
        A
        B
        C
        c = curve;
        c.ext(:,1) = c.ext(:,1)-contact;
        c.ext(:,2) = c.ext(:,2)-zeroF;
        figure
        plot((c.ext(:,1)).*1e6, (real(abs(c.ext(:,2).*1e9).^(1/B)).*abs(c.ext(:,2))./c.ext(:,2)), 'Color', 'red')
        hold on
        plot((ext_lin(:,1)).*1e6, (ext_lin(:,1).*fitobj2.p1 + fitobj2.p2).*(1e9.^(1/B)), 'Color', 'blue', 'linewidth', 1.5)
        xlim([-.1 .5])
        ylim([-.2 .6])
        xlabel('Indentation [um]')
        ylabel('(Force [nN]) ^ 1/B')
        set(gca,'FontSize',12)
        figure
        i0 = get_index(c.ext(:,1), offset);
        x = c.ext(i0:end,1);
        plot((c.ext(:,1)).*1e6,(c.ext(:,2)).*1e9, 'Color', 'red')
        hold on
        plot(x.*1e6, ((x.*fitobj2.p1 + fitobj2.p2).^B).*1e9, 'color', 'blue', 'linewidth', 1.5)
        xlim([-1.2 .5])
        ylim([-.1 .4])
        xlabel('Indentation [um]')
        ylabel('Force [nN]')
        set(gca, 'FontSize', 12)
    end
    
    % Calculate the error from fitting
    sa = gof.a_se;
    sb = gof.b_se;
    sc = gof2.p1_se;
    dF = abs(fc.ext(first,2) - fc.ext(last,2));
    s1(i)=E(i)^2*sa^2/A^2;
    s2(i)=E(i)^2 * log(C)^2*sb^2;
    s3(i) = 0;
    %s3(i)=E(i)^2*log(dF)^2*sb^2/B^2;
    s4(i)=B^2*E(i)^2*sc^2/C^2;
    s(i) = sqrt(s1(i) + s2(i) + s3(i) + s4(i));
    
%     % Optional: LSQ fit for E with fully constrained contact point
%     fn = @(p, d) (1/(1-nu^2)) .* p(1)*1e9 .* lams(d);
%     opts = optimset('Display','off');
%     results = lsqcurvefit(fn, E(i), first:last, fc.ext(first:last,2)'.*1e9, 0, Inf, opts);
%     E(i) = results(1);
%     results(2) = 0;
%     offset(i) = 0;
%     % Optional: LSQ fit for E, semi-constrained contact point in force
%     fn = @(p, d) p(2) + (1/(1-nu^2)) .* p(1)*1e9 .* lams(d);
%     opts = optimset('Display','off');
%     results = lsqcurvefit(fn, [E(i) 0], first:last, fc.ext(first:last,2)'.*1e9, [0 -Inf], [Inf Inf], opts);
%     E(i) = results(1);
%     offset(i) = results(2) * 1e-9;
% %     % Optional: Plot this fit, comment if no LSQ fit
%     figure
%     fn = @(p, d) (1/(1-nu^2)) .* p(1)*1e9 .* lams(d);
%     hold on
%     plot(c.ext(:,1).*1e6, c.ext(:,2).*1e9,'color','red','LineWidth',1)
%     plot(1e6.*[0; fc.ext(:,1)], 1e9.*[offset(i) results(2)/1e9+fn(E(i)./1e9, 1:length(fc.ext))],'color','blue', 'linewidth', 1.5)
%     c = zero_contact(curve,contact);
%     xlabel('Indentation [um]')
%     ylabel('Force [nN]')
%     set(gca, 'FontSize', 16)
    
    % Calculate r2
%     Favg = mean(fc.ext(first:last,2));
%     Fpred = E(i)./(1-nu.^2).*lams(first:last)';
%     Fpred = Fpred - (mean(Fpred)-mean(fc.ext(first:last,2)));
%     sstot = sum((fc.ext(first:last,2) - Favg).^2);
%     ssres = sum((fc.ext(first:last,2) - Fpred).^2);
%     r2(i)=(1-ssres/sstot);
%     % Optional: Plot this fit
%     plot(fc.ext(first:last,1).*1e6, fc.ext(first:last,2).*1e9);
%     hold on
%     plot(fc.ext(first:last,1).*1e6, Fpred.*1e9, 'color', 'red');

end

X = min_ind+binning/2:binning:max_ind-binning/2;
if isempty(X)
    X = (max_ind-min_ind) ./ 2;
end

% Optional: plot the results
if strcmp(doplot,'yes')
%     figure(1)
%     x = min_ind+binning/2:binning:max_ind-binning/2;
%     if isempty(x); x = 0; end
%     E = E(1:end) .* 1e-3;
%     x = x(1:end) .* 1e6;
%     s = s(1:end) .* 1e-3;
%     % errorbar(x,E,s,'LineWidth',1);
%     plot(x,E,'LineWidth',3);
%     set(gca, 'FontSize', 16)
%     xlabel('Indentation (\mum)')
%     ylabel('Apparent Youngs Modulus (kPa)')
%     %hold all
%     E = E * 1e3;
%     fc.ext(:,1) = fc.ext(:,1) * 1e-6;
%     fc.ext(:,2) = fc.ext(:,2) * 1e-9;
    
    % Optional: plot results with force curve in plotyy style
%     x = min_ind+binning/2:binning:max_ind-binning/2;
%     E = real(E);
%     for i=length(x):-1:1
%         if E(i) < 0
%             E(i) = [];
%             x(i) = [];
%             s(i) = [];
%         end
%     end
%     E = E(1:end) .* 1e-3;
%     x = x(1:end) .* 1e6;
%     s = s(1:end) .* 1e-3;
%     fc.ext(:,1) = fc.ext(:,1) * 1e6;
%     fc.ext(:,2) = fc.ext(:,2) * 1e9;
%     [ax,h1,h2] = plotyy(fc.ext(:,1),fc.ext(:,2),x,E, ...
%                 @(x,y)plot(x,y,'LineWidth',3,'Color','red'), ...
%                 @(x,E)plot(x,E,'LineWidth',3,'Color','blue'));
%                 %@(x,E)errorbar(x,E,s,'LineWidth',3,'Color','blue'));
%     axis(ax(1),[eps max(fc.ext(end,1)) eps max(fc.ext(:,2))])
%     %axis(ax(1),[eps max(fc.ext(end,1)) eps 100])
%     % axis(ax(1),[0 12 0 50])
%     axis(ax(2),[eps max(fc.ext(end,1)) min(E)-max(s) max(E)+max(s)])
%     % axis(ax(2),[eps max(fc.ext(end,1)) 0 20])
%     % axis(ax(2),[0 12 .4 1.6])
%     set(get(ax(1),'XLabel'), 'String', 'Indentation (\mum)', 'FontSize', 14)
%     set(get(ax(1),'YLabel'), 'String', 'Force (nN)', 'FontSize', 14)
%     set(get(ax(2),'YLabel'), 'String', 'Apparent Youngs Modulus (kPa)', 'FontSize', 14)
%     set(ax(1), 'YTickMode', 'auto', 'FontSize', 12, 'box','off')
%     set(ax(2), 'YTickMode', 'auto', 'FontSize', 12)
%     hold all
%     E = E * 1e3;
%     fc.ext(:,1) = fc.ext(:,1) * 1e-6;
%     fc.ext(:,2) = fc.ext(:,2) * 1e-9;
    
    % Optional: reconstruct the piecewise force curve and plot it and get sse
%     figure(2)
%     ext_predicted = [];
%     for i=1:nbins
%         first = get_index(fc.ext(:,1), min_ind+binning*(i-1));
%         last = get_index(fc.ext(:,1), min_ind+binning*i);
%         ext = fc.ext(first:last,:);
%         ext_pred = [ext(:,1) (E(i)/(1-nu^2)).*lams(first:last)'];
%         ext_pred(:,2) = ext_pred(:,2) - (mean(ext_pred(:,2)-ext(:,2)));
%         ext_predicted = [ext_predicted; ext_pred];
%         sse(i) = sum(abs(ext_pred(:,2)-ext(:,2)))/length(ext(:,1));
%     end
%     fc2 = zero_contact(curve,contact);
%     fc2 = curve;
%     fc2.ext(:,1) = fc2.ext(:,1) - contact;
%     fc2.ext(:,2) = fc2.ext(:,2) - zeroF;
%     plot(fc2.ext(:,1)*1e6, fc2.ext(:,2)*1e9, 'LineWidth', 3, 'color', 'red')
%     hold on
%     plot(ext_predicted(:,1)*1e6, (ext_predicted(:,2))*1e9, 'LineWidth', 6, 'color', [0 .5 0]);
%     set(gca, 'FontSize', 12)
%     xlabel('Indentation [\mum]')
%     ylabel('Force [nN]')
    
    % Optional: plot the a's to see contact radius as a function of
    % indentation depth
    % plot(fc.ext(:,1)*1e6, as'*1e6, 'LineWidth', 4);
    % xlabel('Indentation [\mum]', 'FontSize', 14);
    % ylabel('Contact Radius [\mum]', 'FontSize', 14);
    % set(gca, 'FontSize', 12)
end

end

