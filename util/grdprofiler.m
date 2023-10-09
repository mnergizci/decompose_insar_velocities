function [pts]=grdprofiler(grdfile, prof_width, bin_size, pts)
%% Function to profile either .grd or .geo.tif files

% Jack McGrath, Aug 2023

if nargin < 2 || prof_width == 0
    prof_width = 50; %km
end

if nargin < 3 || bin_size == 0
    bin_size = 5; %km
end

if exist([grdfile, '.grd'], 'file')
    [lon,lat,data] = grdread2([grdfile, '.grd']);
    fprintf('Loaded %s.grd\n', grdfile)
elseif exist([grdfile, '.geo.tif'], 'file')
    [lon,lat,data,~,~] = read_geotiff([grdfile, '.geo.tif'],'single');
    fprintf('Loaded %s.geo.tif\n', grdfile)
else
    print('Cant find %s as .grd or .geo.tif. Aborting....\n')
    return
end

% Select Profile Location

vik = load('vik.mat');
vik = vik.vik;
figure
imagesc(lon,lat,data,'AlphaData',~isnan(data))
if strcmpi(grdfile(end-1:end), 'vU')
    hold on
    vgps=readmatrix('/nfs/a285/homes/eejdm/IanVels/NZ_InSAR_GPS/NZ_Vertical_GPS_2003-2011_GRL.txt');
    in=inpolygon(vgps(:,1),vgps(:,2),[167,167,174.5,174.5],[-47 -40 -40 -48]);
    vgps=vgps(in,:);
    scatter(vgps(:,1), vgps(:,2), 15, vgps(:,3),'filled', 'MarkerEdgeColor','black')
    hold off
end
axis xy
colorbar
caxis([-prctile(abs(data(:)),95), prctile(abs(data(:)),95)])
pbaspect([1,1,1])
colormap(vik)

if nargin < 4
    fprintf('Select Profile\n')
    pts=ginput;
    fprintf('[%.3f, %.3f; %.3f, %.3f]\n', pts')
else
    fprintf('Using input points\n')
end
hold on
plot(pts(1,1),pts(1,2),'g*')
plot(pts(2,1),pts(2,2),'r*')
plot(pts(:,1), pts(:,2),'k')

% Convert points to UTM

[XX,YY] = meshgrid(lon,lat);
nonan_pix = ~isnan(data);

data=data(nonan_pix);
XX=XX(nonan_pix);
YY=YY(nonan_pix);

llh = [XX,YY];

xy=llh2local(pts',pts(1,:))';
utm=llh2local(llh',pts(1,:))';

xy(:,3) = 0;
utm(:,3) = 0;

half_width = prof_width/2;

% Select points in profile swath
prof_len = sqrt(sum((xy(2,:) - xy(1,:)).^2));
bearing=atand((xy(2,1) - xy(1,1))/ (xy(2,2) - xy(1,2)));

rotmat = rotz(bearing);
xyrot = xy / rotmat;
utmrot = utm / rotmat;
dists = utmrot(:,2) * (abs(xyrot(2,2))/xyrot(2,2));
prof_pix = find(abs(utmrot(:,1)) <= half_width & dists >= 0 & dists < prof_len);

if length(prof_pix) == 0
    boxrot=[-half_width, half_width, half_width, -half_width;0,0,prof_len,prof_len;0,0,0,0]';
    boxrot(:,2) = boxrot(:,2) * (abs(xyrot(2,2))/xyrot(2,2));
    boxutm = boxrot * rotmat;
    profbox=local2llh(boxutm(:,1:2)',pts(1,:))';
    plot(profbox([1,2,3,4,1],1),profbox([1,2,3,4,1],2), 'k')
    
    
    fprintf('No data in profile. Doing Whataroa transect instead\n')
    pts = [169.633, -43.442; 171.168, -44.583];
    fprintf('[%.3f, %.3f; %.3f, %.3f]\n', pts')
    plot(pts(1,1),pts(1,2),'g*')
    plot(pts(2,1),pts(2,2),'r*')
    plot(pts(:,1), pts(:,2),'k')
    
    xy=llh2local(pts',pts(1,:))';
    utm=llh2local(llh',pts(1,:))';
    
    xy(:,3) = 0;
    utm(:,3) = 0;
    
    half_width = prof_width/2;
    
    % Select points in profile swath
    prof_len = sqrt(sum((xy(2,:) - xy(1,:)).^2));
    bearing=atand((xy(2,1) - xy(1,1))/ (xy(2,2) - xy(1,2)));
    
    rotmat = rotz(bearing);
    xyrot = xy / rotmat;
    utmrot = utm / rotmat;
    dists = utmrot(:,2) * (abs(xyrot(2,2))/xyrot(2,2));
    prof_pix = find(abs(utmrot(:,1)) <= half_width & dists >= 0 & dists < prof_len);
end

boxrot=[-half_width, half_width, half_width, -half_width;0,0,prof_len,prof_len;0,0,0,0]';
boxrot(:,2) = boxrot(:,2) * (abs(xyrot(2,2))/xyrot(2,2));
boxutm = boxrot * rotmat;
profbox=local2llh(boxutm(:,1:2)',pts(1,:))';
plot(profbox([1,2,3,4,1],1),profbox([1,2,3,4,1],2), 'k')

data = data(prof_pix);

dists = dists(prof_pix);

if strcmpi(grdfile(end-1:end), 'vU')
    utmgps=llh2local(vgps(:,[1,2])',pts(1,:))';
    utmgps(:,3)=0;
    gpsrot = utmgps / rotmat;
    gpsdists = gpsrot(:,2) * (abs(gpsrot(2,2))/gpsrot(2,2));
    gps_pix = find(abs(gpsrot(:,1)) <= half_width & gpsdists >= 0 & gpsdists < prof_len);
    vgps=vgps(gps_pix,:);
    gpsdists=gpsdists(gps_pix);
elseif strcmpi(grdfile(end-1:end), 'vE')
    gps=readmatrix('../GNSS/beavan_2016_gps_data_synth_SI.csv');
    in=inpolygon(gps(:,1),gps(:,2),[167,167,174.5,174.5],[-47 -40 -40 -48]);
    gps=gps(in,:);
    utmgps=llh2local(gps(:,[1,2])',pts(1,:))';
    utmgps(:,3)=0;
    gpsrot = utmgps / rotmat;
    gpsdists = gpsrot(:,2) * (abs(gpsrot(2,2))/gpsrot(2,2));
    gps_pix = find(abs(gpsrot(:,1)) <= half_width & gpsdists >= 0 & gpsdists < prof_len);
    gps=gps(gps_pix,:);
    gpsdists=gpsdists(gps_pix);
end

% Bin profile data

bin_edges=0:bin_size:prof_len+bin_size;
bin_centers=(bin_size/2):bin_size:prof_len;
n_bins=length(bin_centers);
bin_vals=nan(1, n_bins);
bin_stds=nan(1, n_bins);

for ii=1:n_bins
    in_bin=find(bin_edges(ii) <= dists & dists < bin_edges(ii+1));
    bin_vals(ii) = median(data(in_bin));
    bin_stds(ii) = std(data(in_bin));
end

cmin=floor(prctile(data,1) / 5) * 5;
cmax=ceil(prctile(data,99) / 5) * 5;

figure()
plot(dists, data, '.', 'markersize',2)
hold on
plot(bin_centers, bin_vals, 'k', 'linewidth', 2)
plot(bin_centers, bin_vals + bin_stds, 'k', 'linewidth', 1)
plot(bin_centers, bin_vals - bin_stds, 'k', 'linewidth', 1)
if strcmpi(grdfile(end-1:end), 'vU')
    errorbar(gpsdists, vgps(:,3), vgps(:,4),'r.')
    cmax = ceil(max([cmax,max(vgps(:,3))]) / 5) * 5;
elseif strcmpi(grdfile(end-1:end), 'vE')
    errorbar(gpsdists, gps(:,3), gps(:,5),'r.')
    cmax = ceil(max([cmax,max(gps(:,3))]) / 5) * 5;
end

plot(bin_centers, bin_vals, 'k.', 'markersize', 4)
title(grdfile)

if exist(['hgt.grd'], 'file')
    [hlon,hlat,hgt] = grdread2('hgt.grd');
    fprintf('Loaded %s.grd\n', 'hgt')
    nonan_pix = ~isnan(hgt);
    [XX,YY] = meshgrid(hlon,hlat);
    hgt=hgt(nonan_pix);
    XX=XX(nonan_pix);
    YY=YY(nonan_pix);
    llh = [XX,YY];
    utm=llh2local(llh',pts(1,:))';
    utm(:,3) = 0;
    utmrot = utm / rotmat;
    hgtdists = utmrot(:,2) * (abs(xyrot(2,2))/xyrot(2,2));
    prof_pix = find(abs(utmrot(:,1)) <= half_width & hgtdists >= 0 & hgtdists < prof_len);
    hgtdists=hgtdists(prof_pix);
    hgt=hgt(prof_pix);
    hgt_vals=nan(1, n_bins);
    hgt_stds=nan(1, n_bins);
    for ii=1:n_bins
        in_bin=find(bin_edges(ii) <= hgtdists & hgtdists < bin_edges(ii+1));
        hgt_vals(ii) = nanmedian(hgt(in_bin));
        hgt_stds(ii) = nanstd(hgt(in_bin));
    end
    hmax=(ceil((max(hgt_vals + hgt_stds) / 500)) * 500) * 2;
    cmin=cmin - cmax + cmin;
end

ylim([cmin cmax])
if exist(['hgt.grd'], 'file')
    yyaxis right
    plot(bin_centers, hgt_vals, 'k', 'linewidth', 2)
    plot(bin_centers, hgt_vals + hgt_stds, 'k', 'linewidth', 1)
    plot(bin_centers, hgt_vals - hgt_stds, 'k', 'linewidth', 1)
    ylim([0 hmax])
end
xlabel('Distance (km)')

%% Write outputs
save_flg = input('Save these results? (y/n)\n', "s");

if strcmpi(save_flg, 'y')
    if exist(['hgt.grd'], 'file')
        bindata=[bin_centers; bin_vals; bin_stds; bin_vals-bin_stds; bin_vals+bin_stds; hgt_vals; hgt_stds; hgt_vals-hgt_stds; hgt_vals+hgt_stds];
        pointdata=[dists, data, hgt]';
    else
        bindata=[bin_centers; bin_vals; bin_stds; bin_vals-bin_stds; bin_vals+bin_stds];
        pointdata=[dists, data]';
    end
    
    outdir=sprintf('prof_%s_%.1f_%.1f', grdfile, pts(1,1), pts(2,1));
    mkdir(outdir)
    outfile=[outdir, filesep, grdfile,'_profLine.txt'];
    fid=fopen(outfile, 'w');
    fprintf(fid, '>>> Prof Line %.2f x %.2f km Min: %.2f Max: %.2f\n', prof_len, prof_width, min(data), max(data));
    fprintf(fid, '%.2f %.2f\n',pts');
    fclose(fid);
    
    outfile=[outdir, filesep, grdfile,'_profBox.txt'];
    fid=fopen(outfile, 'w');
    fprintf(fid, '>>> Prof Box\n');
    fprintf(fid, '%.2f %.2f\n',profbox');
    fclose(fid);
    
    outfile=[outdir, filesep, grdfile,'_profBins.txt'];
    fid=fopen(outfile, 'w');
    if exist(['hgt.grd'], 'file')
        fprintf(fid, '>>> Dist; Value; std; val-std; val+std; hgt; hgtstd; hgtval-std; hgtval+std\n');
        fprintf(fid, '%.2f %.2f %.2f %.2f %.2f\n',bindata);
    else
        fprintf(fid, '>>> Dist; Value; std; val-std; val+std\n');
        fprintf(fid, '%.2f %.2f %.2f %.2f %.2f\n',bindata);
    end
    fclose(fid);
    
    outfile=[outdir, filesep, grdfile,'_profPoints.txt'];
    fid=fopen(outfile, 'w');
    if exist(['hgt.grd'], 'file')
        fprintf(fid, '>>> Dist; Value; hgt\n');
        fprintf(fid, '%.2f %.2f %.2f\n',pointdata);
    else
        fprintf(fid, '>>> Dist; Value\n');
        fprintf(fid, '%.2f %.2f\n',pointdata);
    end
    fclose(fid);
    
    if length(gpsdists) > 0
        if strcmpi(grdfile(end-1:end), 'vU')
            GNSS=[gpsdists, vgps(:,3), vgps(:,4)];
        elseif strcmpi(grdfile(end-1:end), 'vE')
            GNSS=[gpsdists, gps(:,3), gps(:,5)];
        end
        
        outfile=[outdir, filesep, grdfile,'_GNSS.txt'];
        fid=fopen(outfile, 'w');
        fprintf(fid, '>>> Dist; Vel; std\n');
        fprintf(fid, '%.2f %.2f %.2f\n',GNSS');
        fclose(fid);
    end
    end