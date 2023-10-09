function [vel, gnss_los] = ref_to_gnss_fields(par,cpt,xx,yy,vel,compE,compN,compU, gnss_E,gnss_N, gnss_U, asc_frames_ind,desc_frames_ind,insar_id, hgt)
%=================================================================
% function ref_to_gnss()
%-----------------------------------------------------------------
% Tie InSAR velocities into a GNSS reference frame.
% Largely based on the method from Weiss et al. (2020).
% Main steps are:
%   - interpolate GNSS stations velocities into continuous fields
%       (provided as an input)
%   - project GNSS into LOS for each frame/track
%   - calculate residual between InSAR and GNSS
%   - either fit a polynomial to the residual, or apply a filter
%   - subtract this smoothed difference from the InSAR
%
% INPUT:
%   par: parameter structure from readparfile.
%   cpt: structure containing colour palettes
%   x, y: vectors of longitude and latitude
%   vel: regridded velocities (3D array)
%   compE, compN, compU: regridded component vectors (3D arrays)
%   vstd: regridded velocity uncertainties
%   asc_frames_ind, desc_frames_ind: indices for ascending and descending
%       frames/tracks
% OUTPUT:
%   vel: velocities in GNSS reference system
% gnss_los: GNSS velocities in LOS for each frame
%
%
% Andrew Watson     06-06-2022
% Jack McGrath      05-07-2023 Output GNSS LOS
% Jack McGrath      26-07-2023 Store and Load Resid Planes
%=================================================================
nframes = size(vel,3);

if ~ismissing(par.use_stored_ref_planes)
    if par.merge_tracks_along == 2
        residfile = [par.out_path, 'GNSS', filesep, 'GNSS_track_resid', par.use_stored_ref_planes, '.mat'];
    else
        residfile = [par.out_path, 'GNSS', filesep, 'GNSS_frame_resid', par.use_stored_ref_planes, '.mat'];
    end
    if isfile(residfile)
        fprintf('Loading GNSS Resid file %s\n', residfile)
        resids = load(residfile);
        
        % Check that everything matches
        if length(resids.insar_id) ~= length(insar_id)
            warning('Different number of InSAR IDs in loaded and processing datasets. Calculating reference from scratch\n')
            par.use_stored_ref_planes = missing;
        elseif sum(ismember(resids.insar_id, insar_id)) ~= length(resids.insar_id)
            warning('Different InSAR IDs in loaded and processing datasets. Calculating reference from scratch\n')
            par.use_stored_ref_planes = missing;
        else
            % InSAR IDs are the same, ensure that they are in the same order
            [~,~,sort_order] = intersect(resids.insar_id,insar_id,'stable');
            if length(sort_order) ~= length(insar_id);
                fprintf('CAUTION: Some along-track frames have not been merged so order cannot be assigned\n')
                fprintf('         Just going to assume that you havent changed the input order instead....\n')
                sort_order = 1:nframes;
            end                
        end
        gnss_resid_plane = resids.gnss_resid_plane(:, :, sort_order);
        gnss_los = resids.gnss_los(:, :, sort_order);
        
    else
        fprintf('Residual file does not exist. Calculating reference from scratch\n')
        par.use_stored_ref_planes = missing;
    end
end

if ismissing(par.use_stored_ref_planes)
    % pre-allocate
        gnss_resid_plane = zeros([size(xx) nframes]);
        gnss_los = zeros([size(xx) nframes]);
end

% coords
x = xx(1,:); y = yy(:,1);

for ii = 1:nframes
    % skip loop if vel is empty (likely because of masking)
    if all(isnan(vel(:,:,ii)),'all')
        disp(['Layer ' num2str(ii) ' of vel is empty after masking, skipping referencing'])
        continue
    end
    if ismissing(par.use_stored_ref_planes)

        % for plotting
        if par.plt_ref_gnss_indv == 1
            vel_orig = vel(:,:,ii);
        end

        % convert gnss fields to los
        if isempty(gnss_U)
            gnss_los(:,:,ii) = (gnss_E.*compE(:,:,ii)) + (gnss_N.*compN(:,:,ii));
        else
            gnss_los(:,:,ii) = (gnss_E.*compE(:,:,ii)) + (gnss_N.*compN(:,:,ii)) + (gnss_U.*compU(:,:,ii));
        end

        % calculate residual
        vel_tmp = vel(:,:,ii);

        % mask after deramping (deramp isn't carried forward)
        % use a hardcoded 10 mm/yr limit to remove large signals (mainly
        % subsidence and seismic)
        vel_deramp = deramp(x,y,vel_tmp);
        vel_deramp = vel_deramp - mean(vel_deramp(:),'omitnan');
        vel_mask = vel_deramp > par.refmask_max | vel_deramp < par.refmask_min;

        vel_tmp(vel_mask) = nan;

        gnss_resid = vel_tmp - gnss_los(:,:,ii);

        if par.remove_linear_APS == 1
            % Find APS component from LOS
            los_resid = gnss_resid;
            gnss_xx = xx(~isnan(los_resid));
            gnss_yy = yy(~isnan(los_resid));
            hgt_tmp = hgt(:,:,ii);
            hgt_tmp = hgt_tmp(~isnan(los_resid));
            los_resid = los_resid(~isnan(los_resid));

            % Calculate linear APS from residuals
            G_aps = [ones(length(gnss_xx),1) hgt_tmp];
            m_aps = (G_aps'*G_aps)^-1*G_aps'*los_resid;
            aps_corr = m_aps(2) .* hgt(:,:,ii);
%             if par.plt_ref_gnss_indv == 1
%             figure(), hold on
%             plot(hgt_tmp,los_resid,'.')
%             plot(hgt_tmp, los_resid - m_aps(2) * hgt_tmp,'.')
%             legend([{'Original'},{'APS Corrected'}])
%             end

            % Subtract this correction from original velocity, reclaculate
            % residual
            vel(:,:,ii) = vel(:,:,ii) - aps_corr;
            vel_tmp = vel(:,:,ii);
            m_aps = m_aps(2);

            % mask after deramping (deramp isn't carried forward)
            % use a hardcoded 10 mm/yr limit to remove large signals (mainly
            % subsidence and seismic)
            vel_deramp = deramp(x,y,vel_tmp);
            vel_deramp = vel_deramp - mean(vel_deramp(:),'omitnan');
            vel_mask = vel_deramp > par.refmask_max | vel_deramp < par.refmask_min;

            vel_tmp(vel_mask) = nan;
            gnss_resid = vel_tmp - gnss_los(:,:,ii);
        end
        
        % method switch
        switch par.ref_type
            case 1 % polynomial surface
                
                % remove nans
                gnss_xx = xx(~isnan(gnss_resid));
                gnss_yy = yy(~isnan(gnss_resid));
                hgt_tmp = hgt(:,:,ii);
                hgt_tmp = hgt_tmp(~isnan(gnss_resid));
                gnss_resid = gnss_resid(~isnan(gnss_resid));
                
                % centre coords
                midx = (max(gnss_xx) + min(gnss_xx))/2;
                midy = (max(gnss_yy) + min(gnss_yy))/2;
                gnss_xx = gnss_xx - midx ;gnss_yy = gnss_yy - midy;
                all_xx = xx - midx; all_yy = yy - midy;
                
                % check that an order has been set
                if isempty(par.ref_poly_order)
                    error('Must set par.ref_poly_order if using poly for referencing')
                end
                
                % fit polynomial
                if par.ref_poly_order == 1 % 1st order
                    G_resid = [ones(length(gnss_xx),1) gnss_xx gnss_yy];
                    m_resid = (G_resid'*G_resid)^-1*G_resid'*gnss_resid;
                    gnss_resid_plane(:,:,ii) = m_resid(1) + m_resid(2).*all_xx + m_resid(3).*all_yy;
                    
                elseif par.ref_poly_order == 2 % 2nd order
                    G_resid = [ones(length(gnss_xx),1) gnss_xx gnss_yy gnss_xx.*gnss_yy ...
                        gnss_xx.^2 gnss_yy.^2];
                    if par.remove_linear_APS == 2
                        G_resid(:,7) = hgt_tmp;
                    end                    
                    m_resid = (G_resid'*G_resid)^-1*G_resid'*gnss_resid;
                    gnss_resid_plane(:,:,ii) = m_resid(1) + m_resid(2).*all_xx + m_resid(3).*all_yy ...
                        + m_resid(4).*all_xx.*all_yy + m_resid(5).*all_xx.^2 + m_resid(6).*all_yy.^2;
                    if par.remove_linear_APS == 2
                        aps_corr = m_resid(7) .* hgt(:,:,ii);
                        vel(:,:,ii) = vel(:,:,ii) - aps_corr;
                        m_aps = m_resid(7);
                    end
                end
                
            case 2 % filtering
                
                % make sure filter size is odd
                if mod(par.ref_filter_window_size,2) ~= 1
                    error('Filter window size must be an odd number')
                end
                
                % filter gnss residual
                windsize = [par.ref_filter_window_size par.ref_filter_window_size];
                gnss_resid_filtered = ndnanfilter(gnss_resid,'rectwin',windsize);
                
                % reapply nans
                %             gnss_resid_filtered(isnan(gnss_resid)) = nan;
                gnss_resid_filtered(isnan(vel(:,:,ii))) = nan;
                
                % store
                gnss_resid_plane(:,:,ii) = gnss_resid_filtered;
        end
    end

    % mask resid with vel (just for plotting)
    vel_mask = single(~isnan(vel(:,:,ii)));
    vel_mask(vel_mask==0) = nan;
    gnss_resid_plane(:,:,ii) = gnss_resid_plane(:,:,ii) .* vel_mask;
    
    % remove from insar
    if par.decomp_method ~= 4
        vel(:,:,ii) = vel(:,:,ii) - gnss_resid_plane(:,:,ii);
    else
        vel(:,:,ii) = vel(:,:,ii) - gnss_los(:, :, ii) - gnss_resid_plane(:,:,ii);
    end
    
    % optional plotting
    if par.plt_ref_gnss_indv == 1
        
        % limits
        clim = [par.plt_cmin par.plt_cmax];
        x = xx(1,:); y = yy(:,1);
        [~,x_ind,y_ind] = crop_nans(vel(:,:,ii),x,y);
        lonlim = x([x_ind(1) x_ind(end)]); latlim = y([y_ind(1) y_ind(end)]);
        
        f = figure();
        if par.remove_linear_APS == 0
            f.Position([1 3 4]) = [600 1600 600];
            tiledlayout(1,3,'TileSpacing','compact');
        else
            f.Position([1 3 4]) = [600 2000 1200];
            tiledlayout(2,3,'TileSpacing','compact');
        end
        
        nexttile; hold on
        imagesc(x,y,vel_orig,'AlphaData',~isnan(vel_orig)); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Original InSAR vel')
        
        if par.remove_linear_APS > 0
            nexttile; hold on
            imagesc(x,y,aps_corr,'AlphaData',~isnan(vel_orig)); axis xy
            xlim(lonlim); ylim(latlim);
            colorbar; colormap(cpt.vik); caxis([-10 10])
            title(sprintf('APS Component (%.2f mm/yr/km)', m_aps*1000))
            
            nexttile; hold on
            hgt_tmp = hgt(:,:,ii);
            plot(hgt_tmp(:), vel_orig(:), '.')
            plot(hgt_tmp(:), vel_orig(:)-aps_corr(:), '.')
            title(sprintf('APS Correction', m_aps*1000))
            legend([{'Original'},{'APS Corrected'}])
            xlabel('Topography (m)')
            ylabel('Velocity (mm/yr)')
            
            nexttile; hold on
            imagesc(x,y,vel_orig - aps_corr,'AlphaData',~isnan(vel_orig)); axis xy
            xlim(lonlim); ylim(latlim);
            colorbar; colormap(cpt.vik); caxis(clim)
            title(sprintf('APS Corrected', m_aps*1000))
        end
        
        nexttile; hold on
        imagesc(x,y,gnss_resid_plane(:,:,ii),'AlphaData',~isnan(gnss_resid_plane(:,:,ii))); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Referencing residual')
        
        nexttile; hold on
        imagesc(x,y,vel(:,:,ii),'AlphaData',~isnan(vel(:,:,ii))); axis xy
        xlim(lonlim); ylim(latlim);
        colorbar; colormap(cpt.vik); caxis(clim)
        title('Referenced InSAR')
        
    end
    
    % report progress
    disp([num2str(ii) '/' num2str(nframes) ' complete'])
    
end

% Fudge for the event that your ENU comps aren't cropped at the coasts
gnss_los(isnan(vel)) = nan;


if par.plt_ref_gnss_surfaces == 1
    % plot referencing functions
    lonlim = [min(x) max(x)];
    latlim = [min(y) max(y)];
    clim = [par.plt_cmin par.plt_cmax];
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'Referencing surfaces')
    
    % ascending tracks
    t(1) = nexttile; hold on
    plt_data(x,y,gnss_resid_plane(:,:,asc_frames_ind),lonlim,latlim,clim,'Ascending (mm/yr)',[],[])
    colormap(t(1),cpt.vik)
    
    % descending tracks
    t(2) = nexttile; hold on
    plt_data(x,y,gnss_resid_plane(:,:,desc_frames_ind),lonlim,latlim,clim,'Descending (mm/yr)',[],[])
    colormap(t(2),cpt.vik)
end

if par.plt_ref_gnss_los == 1
    % plot referencing functions
    lonlim = [min(x) max(x)];
    latlim = [min(y) max(y)];
    clim = [par.plt_cmin par.plt_cmax];
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'GNSS LOS')
    
    % ascending tracks
    t(1) = nexttile; hold on
    plt_data(x,y,gnss_los(:,:,asc_frames_ind),lonlim,latlim,clim,'Ascending (mm/yr)',[],[])
    colormap(t(1),cpt.vik)
    
    % descending tracks
    t(2) = nexttile; hold on
    plt_data(x,y,gnss_los(:,:,desc_frames_ind),lonlim,latlim,clim,'Descending (mm/yr)',[],[])
    colormap(t(2),cpt.vik)
end

if par.grd_ref_gnss_los == 1
    if ~isfolder([par.out_path, 'GNSS'])
        mkdir([par.out_path, 'GNSS'])
    end
    for ii = 1:length(insar_id)
        LOS = gnss_los(:, : ,ii);
        [y, x] = ind2sub(size(LOS), find(~isnan(LOS)));
        ylims=[floor(min(y)/10) * 10, ceil(max(y)/10) * 10];
        xlims=[floor(min(x)/10) * 10, ceil(max(x)/10) * 10];
        if ylims(1) == 0
            ylims(1) = 1;
        end
        if xlims(1) == 0
            xlims(1) = 1;
        end
        fprintf('%.0f/%.0f Writing GNSS LOS for %s to .grd...\n', ii, length(insar_id), insar_id{ii})
        grdwrite2(xx(1, xlims(1):xlims(2)), yy(ylims(1):ylims(2), 1), LOS(ylims(1):ylims(2),xlims(1):xlims(2)), [par.out_path 'GNSS' filesep insar_id{ii} '_GNSS_LOS.grd'])
    end
    grdwrite2(xx(1, :), yy(:, 1), gnss_E(:,:), [par.out_path 'GNSS' filesep 'GNSS_E.grd'])
    grdwrite2(xx(1, :), yy(:, 1), gnss_N(:,:), [par.out_path 'GNSS' filesep 'GNSS_N.grd'])
end

if ~ismissing(par.store_ref_planes)
    if ~isfolder([par.out_path, 'GNSS'])
        mkdir([par.out_path, 'GNSS'])
    end
    if par.merge_tracks_along == 2
        residfile = [par.out_path, 'GNSS', filesep, 'GNSS_track_resid', par.store_ref_planes, '.mat'];
    else
        residfile = [par.out_path, 'GNSS', filesep, 'GNSS_frame_resid', par.store_ref_planes, '.mat'];
    end
    fprintf('Saved GNSS residuals to %s\n', residfile)
    % TO DO: Crop reference frames to InSAR coverage and save as cell array
    % to save space, or allow a list of the InSAR you need to save + apply to
    ref_type = par.ref_type;
    save(residfile, 'gnss_los', 'gnss_resid_plane', 'insar_id', 'ref_type')
end
end