function plot_vertical_results(file_prefix, downsamp, radius)

if nargin < 2
    downsamp = 1;
end

if nargin < 3
    radius = 0.02; % About 2.5km radius
end

outname=[file_prefix, '_vU_comp.png'];


[Lon,Lat,U] = grdread2([file_prefix, '_vU.grd']);

[X,Y] = meshgrid(Lon, Lat);

nonan=~isnan(U);

X=X(nonan);
Y=Y(nonan);
Unon=U(nonan);

vik=load('vik.mat');
vik=vik.vik;

ian=readmatrix('/nfs/a285/homes/eejdm/IanVels/NZ_InSAR_GPS/NZ_InSAR_2003-2011_GRL.txt');
vgps=readmatrix('/nfs/a285/homes/eejdm/IanVels/NZ_InSAR_GPS/NZ_Vertical_GPS_2003-2011_GRL.txt');

in=inpolygon(vgps(:,1),vgps(:,2),[166.5,166.5,174.5,174.5],[-47 -40 -40 -48]);

vgps=vgps(in,:);

TOPS = inpolygon(vgps(:,1),vgps(:,2),[170.7177,171.6094,172.8744,174.8859,174.7615,171.4021,170.2200],[-42.8265,-43.2347,-43.3571,-42.2959,-40.1327,-40.2347,-42.2755]);
midTOPS = inpolygon(vgps(:,1),vgps(:,2),[171.3119  172.4225  173.1099  173.3849  173.6599  173.7128  173.3849  173.1840  172.8244  172.5917  172.1581  171.7773],[-42.8290  -43.1910  -42.8064  -42.4331  -42.0485  -41.8336  -41.4830  -41.4490  -41.2567  -41.8223  -42.0485  -42.4105]);
notmidTOPS=find(midTOPS == 0);


in=inpolygon(ian(:,1),ian(:,2),[166.5,166.5,174.5,174.5],[-47 -40 -40 -48]);

ian=ian(in,:);


f=figure('Position', [5761 1283 1000 1000]);
[ha, pos] = tight_subplot(2,2,[.05 .08],[.03 .03],[.05 .05]);
axes(ha(1))
if downsamp == 1
imagesc(Lon, Lat, U, 'Alphadata', ~isnan(U))
axis xy
else
    scatter(X(1:downsamp:end),Y(1:downsamp:end),2, Unon(1:downsamp:end), 'filled')
end
caxis([-10, 10])
colormap(vik)
c=colorbar;
xlim([166, 174.5])
ylim([-47 -40])
c.Label.String='mm/yr';
pbaspect([0.88,1,1])
title('Sentinel-1 Vertical Velocities')

hold on

scatter(vgps(:,1), vgps(:,2), 15, vgps(:,3),'filled', 'MarkerEdgeColor','black')
scatter(166.5, -40.5, 15, 5,'filled', 'MarkerEdgeColor','black')
text(166.7, -40.5, 'Vertical GNSS')

hold off

axes(ha(2))
scatter(ian(1:downsamp:end,1),ian(1:downsamp:end,2),2, ian(1:downsamp:end,5), 'filled')
caxis([-10, 10])
colormap(vik)
c=colorbar;
xlim([166, 174.5])
ylim([-47 -40])
c.Label.String='mm/yr';
pbaspect([0.88,1,1])
title('Envisat Vertical Velocities')

hold on

scatter(vgps(:,1), vgps(:,2), 15, vgps(:,3),'filled', 'MarkerEdgeColor','black')

hold off

axes(ha(3))

vgps(:,5)=NaN;

SA=[131, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144];

for ii=1:size(vgps,1)
    xref=X - vgps(ii,1);
    yref=Y - vgps(ii,2);
    dist= sqrt(xref.^2 + yref.^2);
    in = find(dist <= radius);
    vgps(ii,5) = nanmedian(Unon(in));
end

F=scatteredInterpolant(X,Y,Unon, 'natural');
projU=F(ian(:,1),ian(:,2));
diff = projU - ian(:,5);
gpsd=vgps(:,5) - vgps(:,3);

scatter(ian(1:downsamp:end,1),ian(1:downsamp:end,2),2, diff(1:downsamp:end), 'filled')
caxis([-10 10])
colormap(vik)
c=colorbar;
xlim([166, 174.5])
ylim([-47 -40])
c.Label.String='mm/yr';
pbaspect([0.88,1,1])
title('Difference (S1 - Envisat)')
hold on
scatter(vgps(:,1), vgps(:,2), 15, gpsd,'filled', 'MarkerEdgeColor','black')

axes(ha(4))

envimean=round(mean(diff),1);
envistd=round(std(diff),1);

gpsmean=round(nanmean(gpsd),1);
gpsstd=round(nanstd(gpsd),1);

SAd = vgps(SA,5) - vgps(SA,3);
SAmean = round(nanmean(SAd),1);
SAstd = round(nanstd(SAd),1);

histogram(diff, -10:0.5:10)
ylabel('Envisat')
hold on
yyaxis right
histogram(vgps(:,5) - vgps(:,3), -10:0.5:10)
%histogram(vgps(TOPS,5) - vgps(TOPS,3), -10:0.5:10)
histogram(vgps(SA,5) - vgps(SA,3), -10:0.5:10)
ylabel('Vertical GNSS')
ylim([0 25])
text(-10, 23.5, ['Envisat: ', num2str(envimean), '\pm', num2str(envistd), ' mm/yr'])
text(-10, 22.0, ['vGNSS: ', num2str(gpsmean), '\pm', num2str(gpsstd), ' mm/yr'])
text(-10, 20.5, ['SA vGNSS: ', num2str(SAmean), '\pm', num2str(SAstd), ' mm/yr'])
hold off
title('Residuals')

saveas(f,outname)
fprintf('%s saved\n', outname)
