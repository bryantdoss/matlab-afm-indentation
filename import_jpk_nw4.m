function curve=import_jpk_nw4 (filename)

% Import JPK data from Nanowizard 4
% Version date: 2019-02-01

% Very poor loading mechanism... hit cancel to be able to choose a
% directory, then pick directory for QI data (already unzipped)
% Loading QI data in Windows may not work
pathname = '';
if nargin < 1
    [filename, pathname] = uigetfile({'*.jpk-force-map;*.jpk-force;*.jpk-qi-data'});
end
if isdir(filename)
    tmpdir = filename;
elseif filename == 0
    tmpdir = uigetdir; % quick hack: if we get cancel then we can load a directory
else
    filename = strcat(pathname,filename);
    if ispc
        tmpdir = '.\jpktmp2';
        unzip(filename, tmpdir);
    else
        tmpdir = '/tmp/jpktmp2';
        system(['unzip ' strrep(filename, ' ', '\ ') ' -d ' tmpdir]);
    end
end

files = dir(fullfile(tmpdir));

forcemap_flag = 0;
for i=1:size(files)
    if strcmp(files(i).name, 'index')
        forcemap_flag = 1;
    end
end

% The channel storage is a little bit weird, need to convert JPK's
% arbitrary measurements for height and vertical deflection into actual
% physical values using offsets and scaling factors.

% Read in data from the header
fid_header = fopen(fullfile(tmpdir,'shared-data','header.properties'));
header = fread(fid_header, inf, 'uint8=>char')';
fclose(fid_header);

% Get the encoding values for the measuredHeight
% Warning: finding the channel number for measuredHeight is a bit bad
index = findstr(header, 'measuredHeight');
line = header(index-24:index+13);
[~,tmp]=strtok(line,'.');
channel_number = strtok(tmp,'.');
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.encoder.scaling.offset']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
zc1 = str2num(dat(2:end));
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.encoder.scaling.multiplier']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
zc2 = str2num(dat(2:end));

% Get the encoding values for vDeflection
% Warning: finding the channel number for vDeflection is a bit bad
index = findstr(header, 'vDeflection');
line = header(index-24:index+13);
[~,tmp]=strtok(line,'.');
channel_number = strtok(tmp,'.');
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.encoder.scaling.offset']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
dc1 = str2num(dat(2:end));
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.encoder.scaling.multiplier']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
dc2 = str2num(dat(2:end));
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.conversion-set.conversion.distance.scaling.offset']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
dc3 = str2num(dat(2:end));
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.conversion-set.conversion.distance.scaling.multiplier']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
dc4 = str2num(dat(2:end));

% Conversion factors for ddist -> force
line = strtok(header(findstr(header, ['lcd-info.' channel_number '.conversion-set.conversion.force.scaling.multiplier']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
k = str2num(dat(2:end));

% Motor position, need to open the other header file
fid_header2 = fopen(fullfile(tmpdir,'header.properties'));
header2 = fread(fid_header2, inf, 'uint8=>char')';
fclose(fid_header2);
line = strtok(header2(findstr(header2, ['environment.xy-scanner-position-map.xy-scanner.motorstage.start-position.x']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
motor_x = str2num(dat(2:end));
line = strtok(header2(findstr(header2, ['environment.xy-scanner-position-map.xy-scanner.motorstage.start-position.y']):end), char(13));
[dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
motor_y = str2num(dat(2:end));


% Now load in all of the curves
cfolders = dir(fullfile(tmpdir,'index'));
cfolders = cfolders(arrayfun(@(x) x.name(1), cfolders) ~= '.');
if forcemap_flag == 0
    cfolders(1).name = '';
else
    for c=1:size(cfolders,1)
        cfolders(c).name = strcat('index/', cfolders(c).name);
    end
end

for c=1:size(cfolders,1)
    if mod(c,2000) == 0
        disp(['Import: ' num2str(c) '/' num2str(size(cfolders,1))])
    end

    % Spring constant of the probe
    curve(c).k = k;
    curve(c).sensitivity = dc4;
    
    % Get the x,y position of the tip
    fid_header = fopen(fullfile(tmpdir,cfolders(c).name,'header.properties'));
    header = fread(fid_header, inf, 'uint8=>char')';
    fclose(fid_header);
    
    % Type of scanning
    line = strtok(header(findstr(header, 'type'):end), char(13));
    [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
    curve(c).scantype = dat(2:end);
    
    % Type of scanning
    line = strtok(header(findstr(header, 'force-settings.type'):end), char(13));
    [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
    curve(c).scantype2 = dat(2:end);
    
    % x-position of the cantilever
    line = strtok(header(findstr(header, 'position.x'):end), char(13));
    [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
    curve(c).tip_x = str2num(dat(2:end));
    
    % y-position of the cantilever
    line = strtok(header(findstr(header, 'position.y'):end), char(13));
    [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
    curve(c).tip_y = str2num(dat(2:end));
    
    % The number of segments in each curve
    if strcmp(curve(c).scantype, 'quantitative-imaging-series') || strcmp(curve(c).scantype2, 'relative-force-settings')
        line = strtok(header(findstr(header, 'force-segments.count'):end), char(13));        
    elseif strcmp(curve(c).scantype, 'force-scan-series')
        line = strtok(header(findstr(header, 'segments.size'):end), char(13));
    end
    
    [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
    curve(c).nsegments = str2num(dat(2:end));
    
    cantilever_z = [];
    cantilever_d = [];
    cantilever_t = [];
    
    % Loop over each segment in the curve
    for s=0:curve(c).nsegments-1
        % Start position of each segment
        curve(c).start(s+1) = length(cantilever_z)+1;
        
        % Segment duration [seconds]
        if strcmp(curve(c).scantype, 'quantitative-imaging-series')
            term = strcat(['extend.duration']);
        elseif strcmp(curve(c).scantype2, 'relative-force-settings')
            term = strcat(['extend-scan-time']);
        elseif strcmp(curve(c).scantype, 'force-scan-series')
            term = strcat(['segment.' num2str(s) '.duration']);
        end
        line = strtok(header(findstr(header, term):end), char(13));
        [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
        curve(c).duration(s+1) = str2num(dat(2:end));
        
        % Segment style
        if strcmp(curve(c).scantype, 'quantitative-imaging-series') || strcmp(curve(c).scantype2, 'relative-force-settings')
            curve(c).style{1} = 'extend';
            curve(c).style{2} = 'retract';
        elseif strcmp(curve(c).scantype, 'force-scan-series')
            term = strcat(['segment.' num2str(s) '.style']);
            line = strtok(header(findstr(header, term):end), char(13));
            [dat,~] = strtok(line, char(10)); [~,dat] = strtok(dat, '=');
            curve(c).style{s+1} = dat(2:end);
        end
        
        % Import height data
        fid_zs = fopen(fullfile(tmpdir,cfolders(c).name,'segments',num2str(s),'channels','measuredHeight.dat'));
        z = fread(fid_zs,inf,'long','ieee-be');
        fclose(fid_zs);
        numpoints = length(z);
        
        % Import deflection data
        fid_ds = fopen(fullfile(tmpdir,cfolders(c).name,'segments',num2str(s),'channels','vDeflection.dat'));
        d = fread(fid_ds,inf,'long','ieee-be');
        fclose(fid_ds);
        
        % Append the data
        if s > 0; starttime = cantilever_t(end); else starttime = 0; end
        cantilever_z = [cantilever_z; ((-z + zc1) .* zc2)];
        cantilever_d = [cantilever_d; (((d + dc1) .* dc2 + dc3) .* dc4)];
        cantilever_t = [cantilever_t; (starttime+((1:numpoints)'.*curve(c).duration(s+1)./numpoints))];
    
    end
    
    % Add everything for saving
    curve(c).zsens = cantilever_z;
    curve(c).defl = cantilever_d;
    curve(c).time = cantilever_t;
    curve(c).motor_x = motor_x;
    curve(c).motor_y = motor_y;
    
    % Force-indentation segment
    for s=1:curve(c).nsegments
        if strcmp(curve(c).style{s}, 'extend')
            if numel(curve(c).start) == 1
                disp(['Warning: only one segment in curve ' num2str(c)])
                z = cantilever_z;
                d = cantilever_d;
            else
                z = cantilever_z(curve(c).start(s):curve(c).start(s+1));
                d = cantilever_d(curve(c).start(s):curve(c).start(s+1));
            end
            curve(c).ext(:,1) = z-d;
            curve(c).ext(:,2) = d.*k;
        end
    end

    % For QI imaging, remove some things to make the files smaller...
    if strcmp(curve(c).scantype, 'quantitative-imaging-series')
        curve(c).zsens = [];
        curve(c).defl = [];
        curve(c).time = [];
    end
    
end

% Remove the temporary directory
if strcmp(tmpdir, '/tmp/jpktmp2') || strcmp(tmpdir, '.\jpktmp2')
    rmdir(tmpdir, 's');
end

end

