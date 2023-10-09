function plot_east_results(file_prefix)

if nargin < 2
    downsamp = 1;
end

outname=[file_prefix, '_vE_comp.png'];

[Lon,Lat,E] = grdread2([file_prefix, '_vE.grd']);

[X,Y] = meshgrid(Lon, Lat);

nonan=~isnan(E);

X=X(nonan);
Y=Y(nonan);
Enon=E(nonan);

vik=load('vik.mat');
vik=vik.vik;

load('../GNSS/gnss_haines.mat');

f=figure('Position', [5761 1283 1000 1000]);
[ha, pos] = tight_subplot(2,2,[.05 .05],[.03 .03],[.05 .07]);
axes(ha(1))
if downsamp == 1
imagesc(Lon, Lat, E, 'Alphadata', ~isnan(E))
axis xy
else
    scatter(X(1:downsamp:end),Y(1:downsamp:end),5, Enon(1:downsamp:end), 'filled')
end
caxis([-40, 40])
colormap(vik)
c=colorbar;
xlim([166, 174.5])
ylim([-47 -40])
c.Label.String='mm/yr';
pbaspect([0.88,1,1])
title('Sentinel-1 East Velocities')


[XX, YY] = ndgrid(gnss_field.x, gnss_field.y);
F=griddedInterpolant(XX,YY,gnss_field.E');

projE = F(X,Y);

projEgrid = nan(size(E));
projEgrid(nonan) = projE;

axes(ha(2))
%scatter(X(1:downsamp:end),Y(1:downsamp:end),5, projE(1:downsamp:end), 'filled')
imagesc(Lon, Lat, projEgrid, 'Alphadata', ~isnan(E))
axis xy
caxis([-40, 40])
colormap(vik)
c=colorbar;
xlim([166, 174.5])
ylim([-47 -40])
c.Label.String='mm/yr';
pbaspect([0.88,1,1])
title('VDoHS East Velocity')

axes(ha(3))

%diff = Enon - projE;
%difflim = max([ceil(abs(prctile(diff,2.5))/5)*5,ceil(abs(prctile(diff,97.5))/5)*5]);
diff = E-projEgrid;
difflim = max([ceil(abs(prctile(diff(:),2.5))/5)*5,ceil(abs(prctile(diff(:),97.5))/5)*5]);
%scatter(X(1:downsamp:end),Y(1:downsamp:end),5, diff(1:downsamp:end), 'filled')
imagesc(Lon, Lat, E-projEgrid, 'Alphadata', ~isnan(E))
axis xy

caxis([-difflim, difflim])
colormap(vik)
c=colorbar;
xlim([166, 174.5])
ylim([-47 -40])
c.Label.String='mm/yr';
pbaspect([0.88,1,1])
title('Difference (S1 - VDoHS)')

axes(ha(4))
diff=diff(:);

difflim = max([ceil(abs(prctile(diff,1))/5)*5,ceil(abs(prctile(diff,99))/5)*5]);
gpsmean=round(nanmean(diff),1);
gpsstd=round(nanstd(diff),1);
gpsmean=nanmean(diff);
gpsstd=nanstd(diff);

histogram(diff, -difflim:0.5:difflim, 'Normalization', 'probability', 'FaceColor', [0.8500 0.3250 0.0980])
ha(4).YAxisLocation = 'right';
ylabel('Probability')
ymax=0.02;
ylim([0 ymax])
text(-difflim * 0.95, ymax * 0.95, sprintf('VDoHS: %.1f\x00B1%.1f mm/yr', gpsmean, gpsstd))   % ['VDoHS: ', num2str(gpsmean), '\pm', num2str(gpsstd), ' mm/yr']
title('Residuals')

saveas(f,outname)
fprintf('%s saved\n', outname)
