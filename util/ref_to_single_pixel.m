function [vel, gnss_E, gnss_sE, gnss_N, gnss_sN, gnss_U, gnss_sU] = ref_to_single_pixel(par,cpt,xx,yy,vel,gnss_E, gnss_sE, gnss_N, gnss_sN, gnss_U, gnss_sU)
%=================================================================
% function ref_to_single_pixel()
%-----------------------------------------------------------------
% Tie InSAR velocities to a single pixel. All frames must cover the same
% pixel
%                                                               
% INPUT:                                                           
%   par: parameter structure from readparfile.
%   cpt: structure containing colour palettes
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN, compU: regridded component vectors (3D arrays)
%   gnss: structure containing "stations" = matrix of gnss station locations, velocities, and uncertainties
%   asc_frames_ind, desc_frames_ind: indices for ascending and descending
%       frames/tracks
% OUTPUT:    
%   vel: velocities in GNSS reference system
%   
% Jack McGrath     10-09-2023
%                                                                  
%=================================================================

nframes = size(vel,3);

for ii = 1:nframes
    
    % get coords for valid pixels
    xx_check = xx(~isnan(vel(:,:,ii)));
    yy_check = yy(~isnan(vel(:,:,ii)));
    
    % distance between all stations and every valid point
    xx_check = repmat(xx_check,1);
    yy_check = repmat(yy_check,1);
    
    dists = sqrt( (xx_check - par.ref_xmin).^2 + (yy_check - par.ref_ymin).^2 );
    
    % number of stations with at least one vel within radius
    n_stations_within = sum(sum(dists <= par.ref_station_radius) >= 1);
    
    if n_stations_within < 1
        error(['Vel ' num2str(ii) ' contains insufficient gnss stations to perform referening'])
    end
    
end

%% perform referencing
    dists = sqrt( (xx - par.ref_xmin).^2 + (yy - par.ref_ymin).^2 );
    within_radius = find(dists <= par.ref_station_radius);


% reference
for ii = 1:nframes
    
    vel_orig = vel(:,:,ii);
    ref_vel = nanmedian(vel_orig(within_radius));
    
    vel(:,:,ii) = vel(:,:,ii) - ref_vel;
    
  
    % optional plotting
    if par.plt_ref_gnss_indv == 1

        % limits
        clim = [par.plt_cmin par.plt_cmax];
        x = xx(1,:); y = yy(:,1);
        [~,x_ind,y_ind] = crop_nans(vel(:,:,ii),x,y);
        lonlim = x([x_ind(1) x_ind(end)]); latlim = y([y_ind(1) y_ind(end)]);
        
        f = figure();
        f.Position([1 3 4]) = [600 1600 600];
        tiledlayout(1,2,'TileSpacing','compact');

        nexttile; hold on
        imagesc(x,y,vel_orig,'AlphaData',~isnan(vel_orig)); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Original InSAR vel')

        nexttile; hold on
        imagesc(x,y,vel(:,:,ii),'AlphaData',~isnan(vel(:,:,ii))); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        plot(par.ref_xmin, par.ref_ymin,'ko', 'MarkerSize', 10, 'MarkerFaceColor','r')
        title(sprintf('Referenced InSAR: {%.2f} mm/yr', ref_vel))
        
    end
    
    % report progress
    disp([num2str(ii) '/' num2str(nframes) ' complete'])
    
end

if ~isempty(gnss_E)
gnss_E = gnss_E - nanmedian(gnss_E(within_radius));
gnss_N = gnss_N - nanmedian(gnss_N(within_radius));
gnss_sE = gnss_sE - nanmedian(gnss_sE(within_radius));
gnss_sN = gnss_sN - nanmedian(gnss_sN(within_radius));
end
if ~isempty(gnss_U)
    gnss_U = gnss_U - nanmedian(gnss_U(within_radius));
    gnss_sU = gnss_sU - nanmedian(gnss_sU(within_radius));
end
    
% report progress
disp(['GNSS complete'])

