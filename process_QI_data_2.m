function ret = process_QI_data_2 (qi_data, filename1, filename2)

slope = nan(1,length(qi_data));
baseline = nan(1,length(qi_data));
c = nan(1,length(qi_data));
start = nan(1,length(qi_data));
for i=1:length(qi_data)
    c(i) = rotation_minimum(qi_data(i));
    start(i) = qi_data(i).ext(1,1);
    baseline(i) = c(i)-qi_data(i).ext(1,1);
end
cutoff = prctile(baseline,98); % find the 2% longest baselines in the dataset
for i=1:length(qi_data)
    if c(i)-qi_data(i).ext(1,1) > cutoff
        [~, slope(i)] = correct_virtual_deflection(qi_data(i), c(i));
    end
end
slope = nanmedian(slope); % find median slope of these for virtual deflection

e = nan(1,length(qi_data));
offset = nan(1,length(qi_data));
curves = qi_data;
for i=1:length(qi_data)
    curves(i) = correct_virtual_deflection(qi_data(i), c(i), slope); % correct virtual deflection
    [e(i),~,~,offset(i)] = dd_fit(curves(i), c(i), 35e-9, 35/2*pi/180, 0.5, 1, 'parabolic', 0, 1, 'no', 0); % fit initial E
end
offset = real(offset);
E_cutoff = 65535; % set cutoff stiffness for plastic, maximum E is 2^16 Pa
disp('Flatten for the local height')
h = -flatten(curves, c+offset, e, E_cutoff); % flatten the contact points for the height at each pixel

e2 = nan(1,length(qi_data));
offset2 = nan(1,length(qi_data));
r2 = nan(1,length(qi_data));
if isempty(excise(h))
    e2 = e;
    disp('Warning: cannot correct for sample thickness')
else 
    for i=1:length(curves)
        [e2(i),~,r2(i),offset2(i)] = dd_fit(curves(i), c(i)+offset(i), 35e-9, 35/2*pi/180, 0.5, 1, 'parabolic', 0e-9, 800e-9, 'no', h(i)); % re-fit the Young's modulus taking height into account
        if i==1955
            [e2_ex,~,r2_ex,offset_ex]=dd_fit(curves(i), c(i)+offset(i), 35e-9, 35/2*pi/180, 0.5, 1, 'parabolic', 0e-9, 800e-9, 'yes', h(i))
            h(i)
        end
    end
end

% cutoffs for the stiffness at each piel
e2=real(e2);
e2(e>E_cutoff)=nan; % cutoff should be plastic
e2(e2>E_cutoff)=nan; % cutoff should be plastic
e2(e<0)=nan; % get rid of less than zero
e2(h<100e-9)=nan; % very shallow is probably also plastic
% e2((c+offset-start)<10e-9)=nan; % bad contact point / no baseline

disp(['Median E: ' num2str(nanmedian(e2))]);
disp(['Mean E: ' num2str(nanmean(e2))]);

% make plots
x = [curves.tip_x];
xrange = max(x)-min(x);
y = [curves.tip_y];
yrange = max(y)-min(y);
l = sqrt(length(curves)) - 1;
[xq,yq] = meshgrid(min(x):xrange/l:max(x)+xrange/(l*3), min(y):yrange/l:max(y)+yrange/(l*3));

[vx,vy,m] = griddata(x, y, e2, xq, yq);
regions = regionprops(~isnan(m), 'PixelList');
len = 100; % domains smaller than 100 pixels should be removed
for a=1:length(regions)
    if size(regions(a).PixelList,1) < len
        for b=1:size(regions(a).PixelList,1)
            m(regions(a).PixelList(b,2),regions(a).PixelList(b,1)) = nan;
        end
    end
end
m_E = m(:);

% % Plot the stiffness map
figure('units','pixels',...
    'position',[100 100 size(m,2)*2 size(m,1)*2],...
    'resize','off','name','Stiffness Map');
a1 = axes('unit', 'pix', 'position', [1 1 size(m,2)*2 size(m,1)*2]);
h1 = pcolor(vx, vy, m);
shading flat
axis off
axis ij
caxis([0 20e3]) % Caco2
% imwrite(uint16(m), filename1);

% Calculate heights
h_new = h;
% h_new(isnan(e2))=nan;
% h_new = 1:length(h_new);
[vx,vy,m] = griddata(x, y, h_new, xq, yq);
m(isnan(m_E)) = nan;
m_h = m(:);

% % Plot the height map
figure('units','pixels',...
    'position',[100 100 size(m,2)*2 size(m,1)*2],...
    'resize','off','name','Height Map');
a1 = axes('unit', 'pix', 'position', [1 1 size(m,2)*2 size(m,1)*2]);
h1 = pcolor(vx, vy, m);
shading flat
colormap(jet)
axis off
axis ij
caxis([0 6e-6]) % Caco2
% imwrite(uint16(m.*1e9), filename2);

% Set return value
ret = [m_h'; m_E'];

end
