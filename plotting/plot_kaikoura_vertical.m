function plot_vertical_results(file_prefix, downsamp, radius)

if nargin < 2
    downsamp = 1;
end

if nargin < 3
    radius = 0.1;
end

outname=[file_prefix, '_comp.png'];

[Lon,Lat,U] = grdread2([file_prefix, '_vU.grd']);
[Lon,Lat,E] = grdread2([file_prefix, '_vE.grd']);

[X,Y] = meshgrid(Lon, Lat);

nonan=~isnan(U);

X=X(nonan);
Y=Y(nonan);
Unon=U(nonan);
Enon=E(nonan);

vik=load('vik.mat');
vik=vik.vik;

gps=readmatrix('/nfs/a1/eejdm/DIV/GNSS/Kaikoura_GNSS.txt');

invgps = zeros(size(gps,1),3) * NaN;
% figure
% hold on
for ii=1:size(gps,1)
    xref=X - gps(ii,1);
    yref=Y - gps(ii,2);
    dist= sqrt(xref.^2 + yref.^2);
    in = find(dist <= radius);
    invgps(ii,1) = nanmedian(Enon(in));
    invgps(ii,3) = nanmedian(Unon(in));
%     plot(X(in),Y(in),'.')
end
compare_sites = find(~isnan(invgps(:,1)));

f=figure('Position', [5761 1283 1500 1000]);
[ha, pos] = tight_subplot(2,2,[.05 .05],[.03 .03],[.05 .05]);
axes(ha(1))
if downsamp == 1
imagesc(Lon, Lat, E, 'Alphadata', ~isnan(E))
axis xy
else
    scatter(X(1:downsamp:end),Y(1:downsamp:end),5, Enon(1:downsamp:end), 'filled')
end
caxis([-2500, 2500])
colormap(vik)
c=colorbar;
xlim([170.5, 174.5])
ylim([-44 -40])
c.Label.String='mm';
pbaspect([0.88,1,1])
title('East Displacements')

hold on

quiver(gps(:,1),gps(:,2), gps(:,3), gps(:,4),'black')
quiver(gps(:,1),gps(:,2), invgps(:,1), gps(:,4),'red')
legend({'GNSS'; 'Inverted'})
hold off

axes(ha(2))
if downsamp == 1
imagesc(Lon, Lat, U, 'Alphadata', ~isnan(U))
axis xy
else
    scatter(X(1:downsamp:end),Y(1:downsamp:end),5, Unon(1:downsamp:end), 'filled')
end
caxis([-1000, 1000])
colormap(vik)
c=colorbar;
xlim([170.5, 174.5])
ylim([-44 -40])
c.Label.String='mm';
pbaspect([0.88,1,1])
title('Vertical Displacement')

hold on

scatter(gps(:,1), gps(:,2), 15, gps(:,5),'filled', 'MarkerEdgeColor','black')
quiver(gps(:,1),gps(:,2), zeros(size(gps,1),1), gps(:,5),'black')
quiver(gps(:,1),gps(:,2), zeros(size(gps,1),1), invgps(:,3),'red')
hold off



axes(ha(3))
histogram(invgps(compare_sites,1) - gps(compare_sites,3), -5000:50:5000)
ymax=25;
text(-5000, ymax*5/6, ['Difference: ', num2str(round(nanmedian(invgps(:,1) - gps(:,3)),1)), '\pm', num2str(round(nanstd(invgps(compare_sites,1) - gps(compare_sites,3)),1)), ' mm'])
ylim([0 ymax])
hold off
title('East Residuals')

axes(ha(4))
histogram(invgps(compare_sites,3) - gps(compare_sites,5), -1000:10:1000)
ymax=25
text(-1000, ymax*5/6, ['Difference: ', num2str(round(nanmedian(invgps(:,3) - gps(:,5)),1)), '\pm', num2str(round(nanstd(invgps(compare_sites,3) - gps(compare_sites,5)),1)), ' mm'])
ylim([0 ymax])
hold off
title('Vertical Residuals')

saveas(f,outname)
fprintf('%s saved\n', outname)
