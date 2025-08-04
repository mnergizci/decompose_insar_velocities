function decompose_insar_velocities(config_file)
%% decompose_insar_velocities.m
%
% Decomposes line-of-sight velocities (primarily from LiCSBAS) into a
% combination of North, East, and Up velocities, assuming no covariance
% between adjacent pixels. All frames are interpolated onto a new grid.
%
% InSAR velocities may be referenced to a GNSS velocity field (ref2gnss)
% via various methods, and then decomposed into East and Vertical
% components.
%
% Functions are include to correct for plate motion bias, to merge LOS
% velocities along-track, and to perform downsampling.
%
% All options and inputs are read from a config file, an example for which
% can be found in this directory.
%
% File formats can be found in the README and in the headers of each
% script.
%
% Written for matlab/2020a.
%
% Andrew Watson     26-04-2021 (Original version)
%%
close all
clear;clc
addpath(genpath('C:\Users\mmuha\OneDrive - University of Leeds\1.Phd_project\decompose_insar_velocities'));
%%

if nargin < 1
    config_file = ['.', filesep, 'EAF_TRrange.conf'];
end

disp('Beginning run')

% add subdirectory paths
addpath util plotting

% begin timer
tic

%% read parameter file

disp('Loading parameter file')

[par,insarpar] = readparfile(config_file);

%% check that inputs exist

disp('Checking that the requested inputs exist')

if insarpar.ninsarfile == 0
    error('No insar files provided - have you remembered to add "framedir:" before each path?')
end

% pre-allocate
to_remove = false(1,length(insarpar.dir));

for ii = 1:length(insarpar.dir)
    
    % get name of velocity file
    namestruct = dir([insarpar.dir{ii} '*' insarpar.id_vel '*']);
    
    % test existance and record directory if no file exists
    if isempty(namestruct)   
        disp(['Removing ' insarpar.dir{ii}])
        to_remove(ii) = true;        
    end
    
end

% remove missing file dirs
insarpar.dir(to_remove) = []; clear to_remove

%% load inputs

disp('Loading inputs')

% number of velocity maps inputted
nframes = length(insarpar.dir);

if length(unique(insarpar.dir)) ~= nframes
    error('At least one framedir has been repeated, please check your config file.')
end

% load main inputs
[lon,lat,dx,dy,lon_comp,lat_comp,vel,vstd,compE,compN,compU,mask,poly_mask,frames,...
    asc_frames_ind,desc_frames_ind,fault_trace,gnss,borders, hgt] = load_inputs(par,insarpar);

% colour palettes  (https://www.fabiocrameri.ch/colourmaps/)
load('vik.mat');
load('batlow.mat');
cpt.vik = vik; cpt.batlow = batlow; 
clear vik batlow

%% preview inputs

if par.plt_input_vels == 1
    
    disp('Plotting preview of input velocities')
    
    plt_asc_desc_cells(par,lon,lat,vel,mask,asc_frames_ind,...
        desc_frames_ind,cpt.vik,[par.plt_cmin par.plt_cmax],borders,'Input velocities')
    
end

%% downsample unit vectors if required
% ENU from licsbas may not be downsampled to the same level as vel by
% default, so check the sizes and downsample if neccessary. 
disp('Checking if unit vectors need downsampling')

for ii = 1:nframes
    
    % calculate downsample factor from matrix size
    comp_dsfac = round(size(compE{ii},1)./size(vel{ii},1));
    
    if comp_dsfac > 1
        
        disp(['Downsampling unit vectors : ' insarpar.dir{ii}])
    
        % downsample components
        [compE{ii},lon_comp{ii},lat_comp{ii}] = downsample_array(compE{ii},...
            comp_dsfac,comp_dsfac,'mean',lon_comp{ii},lat_comp{ii});
        [compN{ii},~,~] = downsample_array(compN{ii},comp_dsfac,comp_dsfac,'mean');
        [compU{ii},~,~] = downsample_array(compU{ii},comp_dsfac,comp_dsfac,'mean');
        
    end

end

disp('Done')

%% downsample

if par.ds_factor > 0
    
    disp(['Downsampling inputs by a factor of ' num2str(par.ds_factor)])
    
    [lon,lat,vel,vstd,compE,compN,compU,mask] ...
        = downsamp(par,nframes,lon,lat,vel,vstd,compE,compN,compU,mask);
    
end

%% scale velocity uncertainties

if par.scale_vstd == 1
    
    disp(['Scaling velocity uncertainties using ' par.scale_vstd_model ' model'])
    
    scale_vstd_misfit = nan(nframes,1);
    
    for ii = 1:nframes
        
        % temp apply mask if request
%         vstd_temp = vstd{ii};
%         if par.use_mask == 1
%             vstd_temp(mask{ii}==0) = NaN;
%         end            
        
        % apply scaling
        [vstd{ii},scale_vstd_misfit(ii)] = scale_vstd(par,lon{ii},lat{ii},vstd{ii},frames{ii});
        
        % report progress
        if (mod(ii,round(nframes./10))) == 0
            disp([num2str(round((ii./nframes)*100)) '% completed']);
        end
    end
    
    % plot all if requested
    if par.plt_scale_vstd_all == 1       
        disp('Plotting scaled uncertainties')        
        plt_asc_desc_cells(par,lon,lat,vstd,mask,asc_frames_ind,...
            desc_frames_ind,cpt.batlow,[0 3],borders,'Scaled uncertainties')      
        
    end
    
end

%% unify grids

disp('Unifying grids')

[x_regrid,y_regrid,xx_regrid,yy_regrid,vel_regrid,vstd_regrid,...
    mask_regrid,compE_regrid,compU_regrid,compN_regrid,gnss_E,gnss_N,gnss_sE,gnss_sN, gnss_U, gnss_sU, hgt] ...
    = unify_grids(par,lon,lat,lon_comp,lat_comp,dx,dy,vel,vstd,mask,compE,compU,compN,gnss, hgt);

% clear ungridded data
clear vel vstd compE compN compU

disp('Grid unification complete')

%% apply mask

if par.use_mask == 1
    disp('Applying masks')    
    [vel_regrid,vstd_regrid,mask_regrid] ...
        = apply_masks(par,x_regrid,y_regrid,vel_regrid,vstd_regrid,...
        mask_regrid,asc_frames_ind,desc_frames_ind,fault_trace,borders);
    
elseif par.use_mask == 0
    % standarise nan values across arrays (mask does this if used)
    vel_regrid(vel_regrid==0) = nan;
    vstd_regrid(isnan(vel_regrid)) = nan;
    
end

%% merge frames along-track (and across-track)

if par.merge_tracks_along > 0
    
    disp('Merging frames along-track')
    
    [vel_regrid,compE_regrid,compN_regrid,compU_regrid,vstd_regrid,tracks, hgt] ...
        = merge_frames_along_track(par,cpt,x_regrid,y_regrid,vel_regrid,...
        frames,compE_regrid,compN_regrid,compU_regrid,vstd_regrid, hgt);
    
    % update number of frames and indexes if frames have been merged
    if par.merge_tracks_along == 2
        nframes = size(vel_regrid,3);
        asc_frames_ind = find(cellfun(@(x) strncmp('A',x(4),4), tracks));
        desc_frames_ind = find(cellfun(@(x) strncmp('D',x(4),4), tracks));
    end
    
    % run across-track merging
    if par.merge_tracks_across > 0        
        disp('Merging frame across-track. Not used in decomposition.')       
        merge_frames_across_track(par,x_regrid,y_regrid,vel_regrid,tracks,...
            compE_regrid,compU_regrid,vstd_regrid)
    end
    
end

%% reference frame bias correction
% remove "reference frame bias" caused by rigid plate motions.

if par.plate_motion == 1
    disp('Applying plate motion correction')
    [vel_regrid] = plate_motion_bias(par,cpt,x_regrid,y_regrid,vel_regrid,...
        compE_regrid,compN_regrid,asc_frames_ind,desc_frames_ind);
end

%% tie to gnss
% Shift InSAR velocities into the same reference frame as the GNSS
% velocities. Method is given by par.ref2gnss. 

if par.merge_tracks_along == 2
    outdirs = tracks;
else
    outdirs = frames;
end

if par.ref2gnss == 1
    disp('Referencing InSAR to GNSS station velocities')
    vel_regrid = ref_to_gnss_stations(par,cpt,xx_regrid,yy_regrid,vel_regrid,...
        compE_regrid,compN_regrid,gnss,asc_frames_ind,desc_frames_ind,outdirs);   
    
elseif par.ref2gnss == 2
    disp('Referencing InSAR to interpolated GNSS velocities')
    [vel_regrid, gnss] = ref_to_gnss_fields(par,cpt,xx_regrid,yy_regrid,vel_regrid,...
        compE_regrid,compN_regrid,compU_regrid, gnss_E,gnss_N, gnss_U,asc_frames_ind,desc_frames_ind,outdirs, hgt);
    
elseif par.ref2gnss == 3
    disp('Referencing InSAR to common pixel')
    [vel_regrid, gnss_E, gnss_sE, gnss_N, gnss_sN, gnss_U, gnss_sU] = ...
        ref_to_single_pixel(par,cpt,xx_regrid,yy_regrid,vel_regrid,gnss_E, gnss_sE, gnss_N, gnss_sN, gnss_U, gnss_sU);
    
end

%% preview inputs (again)

if par.plt_input_vels == 1 && par.decomp_method == 4
    
    disp('Plotting preview of gnss stripped input velocities')
        plt_asc_desc_cells(par,lon,lat,vel_regrid,mask_regrid,asc_frames_ind,...
        desc_frames_ind,cpt.vik,[par.plt_cmin par.plt_cmax],borders,'GNSS Stripped velocities')
elseif par.plt_input_vels == 1
    disp('Plotting preview of gnss referenced input velocities')
    
    plt_asc_desc_cells(par,lon,lat,vel_regrid,mask,asc_frames_ind,...
        desc_frames_ind,cpt.vik,[par.plt_cmin par.plt_cmax],borders,'GNSS Referenced Velocities',1)
end


save_for_comparison=1;
if save_for_comparison == 1 && length(insarpar.dir) == 1
    fprintf('Saving the referenced InSAR LOS\n')
lonlim = [min(x_regrid) max(x_regrid)];
latlim = [min(y_regrid) max(y_regrid)];
grdwrite2(x_regrid,y_regrid,vel_regrid,[insarpar.dir{1}, filesep, 'LOSref.grd'])
grdwrite2(x_regrid,y_regrid,gnss,[insarpar.dir{1}, filesep, 'gnssLOS.grd'])
grdwrite2(x_regrid,y_regrid,vel_regrid - gnss,[insarpar.dir{1}, filesep, 'LOSref_gnss_diff.grd'])
end
%% calculate frame overlaps if requested

if par.frame_overlaps == 1
    disp('Calculating frame overlap statistics')
    if par.merge_tracks_along == 2
        frame_overlap_stats(vel_regrid,tracks,compU_regrid);
    else
        frame_overlap_stats(vel_regrid,frames,compU_regrid);
    end
end

%% check data coverage

% check for at least one ascending and one descending velocity
asc_coverage = any(~isnan(vel_regrid(:,:,asc_frames_ind)),3);
desc_coverage = any(~isnan(vel_regrid(:,:,desc_frames_ind)),3);
both_coverage = all(cat(3,asc_coverage,desc_coverage),3);

% check North GNSS vels if used
if par.decomp_method ~= 3
    gnss_coverage = ~isnan(gnss_N);
    both_coverage = all(cat(3,both_coverage,gnss_coverage),3);
end

%% velocity decomposition

disp('Performing velocity decomposition')

% check that at least one point has at least two look directions
if sum(both_coverage(:)) == 0
    error('No pixels with multiple look directions - did you provide more than one track?')
end

if par.decomp_method == 0 || par.decomp_method == 1 || par.decomp_method == 3 || par.decomp_method == 4
    
    % either remove gnss N, and decompose into vE and vU, or include gnss N
    % in the linear problem, and solve for vE, vU, and vN
    disp('Solving directly for vE and vU')
    [m_east,m_up,var_east,var_up,model_corr,condG_threshold_mask,var_threshold_mask] ...
        = vel_decomp(par,vel_regrid,vstd_regrid,compE_regrid,compN_regrid,...
        compU_regrid,gnss_N,gnss_sN,gnss_E,gnss_sE, both_coverage);

elseif par.decomp_method == 2 || par.decomp_method == 5
    
    % decompose into vE and vUN, then split vUN into vU and vN (Qi's
    % method)
    disp('Solving for vE and vNU as intermediary step')
    [m_east,m_up,var_east,var_up,model_corr,condG_threshold_mask,var_threshold_mask] ...
        = vel_decomp_vE_vUN(par,vel_regrid,vstd_regrid,compE_regrid,compN_regrid,...
        compU_regrid,gnss_N,gnss_sN,both_coverage);
end

%% plot output velocities

par.out_prefix=sprintf('%s_ref%.0f-%.0f_decomp%.0f',par.out_prefix,par.ref2gnss,par.ref_type,par.decomp_method);

if par.remove_linear_APS == 1
    par.out_prefix = [par.out_prefix,'_APS'];
end

if par.plate_motion == 1
    par.out_prefix = [par.out_prefix,'_PM'];
end

if par.use_mask == 1
    par.out_prefix = [par.out_prefix,'_mask'];
end

disp('Plotting decomposed velocities')

lonlim = [min(x_regrid) max(x_regrid)];
latlim = [min(y_regrid) max(y_regrid)];
clim = [par.plt_cmin par.plt_cmax];
vlim = [par.plt_vmin par.plt_vmax];
p10 = prctile(m_up(~isnan(m_up)), 10);
p90 = prctile(m_up(~isnan(m_up)), 90);
vlim = [p10, p90];

% East and vertical
f = figure();
f.Position([1 3 4]) = [600 2800 1200];
t = tiledlayout(1,2,'TileSpacing','compact');
title(t,'Decomposed velocities')

t(1) = nexttile; hold on
if par.decomp_method == 5
    vtitle = 'vUN component (mm/yr)';
else
    vtitle = 'Vertical (mm/yr)';
end
plt_data(x_regrid,y_regrid,m_up,lonlim,latlim,vlim,vtitle,fault_trace,borders)
colormap(t(1),cpt.vik)

t(2) = nexttile; hold on
p10 = prctile(m_east(~isnan(m_east)), 10);
p90 = prctile(m_east(~isnan(m_east)), 90);
clim = [p10, p90];
plt_data(x_regrid,y_regrid,m_east,lonlim,latlim,clim,'East (mm/yr)',fault_trace,borders)
colormap(t(2),cpt.vik)

if par.save_png == 1
    if ~isfolder(par.out_path)
        mkdir(par.out_path)
    end
    if par.decomp_method == 5
        saveas(f, [par.out_path par.out_prefix '_vEvUN.png'])
    else
        saveas(f, [par.out_path par.out_prefix '_vEvU.png'])
    end
end

% North and  coverage
f = figure();
f.Position([1 3 4]) = [600 2800 1200];
t = tiledlayout(1,2,'TileSpacing','compact');
title(t,'Decomposed velocities and coverage')

t(3) = nexttile; hold on
plt_data(x_regrid,y_regrid,gnss_N,lonlim,latlim,[],'North (mm/yr)',fault_trace,borders)

t(4) = nexttile; hold on
coverage = int8(sum(~isnan(vel_regrid),3));
plt_data(x_regrid,y_regrid,coverage,lonlim,latlim,[],'Data coverage',fault_trace,borders)

if par.save_png == 1
    saveas(f, [par.out_path par.out_prefix '_vNcov.png'])
end

%% plot velocity uncertainties

if par.plt_decomp_uncer == 1
    
    disp('Plotting decomposed velocity uncertainties')

    lonlim = [min(x_regrid) max(x_regrid)];
    latlim = [min(y_regrid) max(y_regrid)];
    clim = [0 ceil(prctile([var_up(:); var_east(:)],95) / 2) * 2];

    f = figure();
    f.Position([1 3 4]) = [600 2000 800];
    t = tiledlayout(1,3,'TileSpacing','compact');
    title(t,'Decomposed velocity uncertainties')

    t(1) = nexttile; hold on
    plt_data(x_regrid,y_regrid,var_up,lonlim,latlim,clim,vtitle,fault_trace,borders)
    colormap(t(1),cpt.batlow)

    t(2) = nexttile; hold on
    plt_data(x_regrid,y_regrid,var_east,lonlim,latlim,clim,'East (mm/yr)',fault_trace,borders)
    colormap(t(2),cpt.batlow)
    
    t(3) = nexttile; hold on
    plt_data(x_regrid,y_regrid,model_corr,lonlim,latlim,[-1 1],'Correlation (mm/yr)',fault_trace,borders)
    colormap(t(3),cpt.batlow)
    
    if par.save_png == 1
        if par.decomp_method == 5
            saveas(f, [par.out_path par.out_prefix '_sEsUN.png'])
        else
            saveas(f, [par.out_path par.out_prefix '_sEsU.png'])
        end
    end

end

%% plot variance and cond(G) threshold masks if used

if par.plt_threshold_masks == 1
    
    disp('Plotting cond(G) and variance masks')
    
    f = figure();
    f.Position([1 3 4]) = [600 1600 600];
    
    
    t = tiledlayout(1,2,'TileSpacing','compact');
    title(t,'Threshold masks')
    
    t(1) = nexttile; hold on
    plt_data(x_regrid,y_regrid,condG_threshold_mask,lonlim,latlim,[],'cond(G) mask',fault_trace,borders)
    
    t(2) = nexttile; hold on
    plt_data(x_regrid,y_regrid,var_threshold_mask,lonlim,latlim,[],'variance mask',fault_trace,borders)
    
end

%% save outputs

% save geotiffs
if par.save_geotif == 1
    
    % create out directory
    if ~isfolder(par.out_path)
        mkdir(par.out_path)
    end
    
    disp('Saving the following outputs:')
    disp([par.out_path par.out_prefix '_vE.geo.tif'])
    if par.decomp_method == 5
        disp([par.out_path par.out_prefix '_vUN.geo.tif'])
    else
        disp([par.out_path par.out_prefix '_vU.geo.tif'])
        disp([par.out_path par.out_prefix '_vN.geo.tif'])
    end
    
    % create georeference
    georef = georefpostings([min(y_regrid) max(y_regrid)],...
        [min(x_regrid) max(x_regrid)],size(m_up),'ColumnsStartFrom','south',...
        'RowsStartFrom','west');
    
    % write geotifs
    
    geotiffwrite([par.out_path par.out_prefix '_vE.geo.tif'],m_east,georef)
    geotiffwrite([par.out_path par.out_prefix '_sE.geo.tif'],var_east,georef)
    
    if par.decomp_method == 5
        geotiffwrite([par.out_path par.out_prefix '_vUN.geo.tif'],m_up,georef)
        geotiffwrite([par.out_path par.out_prefix '_sUN.geo.tif'],var_up,georef)
    else
        geotiffwrite([par.out_path par.out_prefix '_vU.geo.tif'],m_up,georef)
        geotiffwrite([par.out_path par.out_prefix '_sU.geo.tif'],var_up,georef)
    end
    
    if exist('m_north','var')
        geotiffwrite([par.out_path par.out_prefix '_vN.geo.tif'],m_north,georef)
    end
    
    copyfile(config_file, [par.out_path par.out_prefix '.conf'])
end

% save grd files
if par.save_grd == 1
    
    % create out directory
    if ~isfolder(par.out_path)
        mkdir(par.out_path)
    end
    
    disp('Saving the following outputs:')
    disp([par.out_path par.out_prefix '_vE.grd'])
    if par.decomp_method == 5
        disp([par.out_path par.out_prefix '_vUN.grd'])
    else
        disp([par.out_path par.out_prefix '_vU.grd'])
    end
    disp([par.out_path par.out_prefix '_mask.grd'])
    
    % write grd files
    grdwrite2(x_regrid,y_regrid,m_east,[par.out_path par.out_prefix '_vE.grd']);
    grdwrite2(x_regrid,y_regrid,var_east,[par.out_path par.out_prefix '_sE.grd']);
    %grdwrite2(x_regrid,y_regrid,mask,[par.out_path par.out_prefix '_mask.grd']);
    
    if par.decomp_method == 5
        grdwrite2(x_regrid,y_regrid,m_up,[par.out_path par.out_prefix '_vUN.grd'])
        grdwrite2(x_regrid,y_regrid,var_up,[par.out_path par.out_prefix '_sUN.grd'])
    else
        grdwrite2(x_regrid,y_regrid,m_up,[par.out_path par.out_prefix '_vU.grd'])
        grdwrite2(x_regrid,y_regrid,var_up,[par.out_path par.out_prefix '_sU.grd'])
    end
    
    if exist('m_north','var')
        disp([par.out_path par.out_prefix '_vN.grd'])
        grdwrite2(x_regrid,y_regrid,m_north,[par.out_path par.out_prefix '_vN.grd']);
    end
    
end
%% save frames
% save frames / tracks
if par.save_frames == 1
    % toggle between tracks and frames depending on if the merge has
    % happened
    if par.merge_tracks_along == 2
        outdirs = tracks;
        
%         % identify repeated tracks and append a number
%         [~,~,rep_ind] = unique(outdirs);
%         for ii = 1:length(tracks)
%             track_ind = strcmp(tracks,tracks{ii});
    
        % THIS NEEDS FINISHING - ERROR WHEN TRACKS GET SPLIT.
            
    else
        outdirs = frames;
    end
    
    % loop through frames / tracks
    for ii = 1:nframes
        
        % create frame directory
        if ~isfolder([par.out_path outdirs{ii}])
            mkdir([par.out_path outdirs{ii}])
        end

        % crop to minimise file size
        [~,x_ind,y_ind,x_crop,y_crop] = crop_nans(vel_regrid(:,:,ii),x_regrid,y_regrid);
        
        % create georeference
        georef = georefpostings([min(y_crop) max(y_crop)],...
            [min(x_crop) max(x_crop)],size(vel_regrid(y_ind,x_ind,:)),'ColumnsStartFrom','south',...
            'RowsStartFrom','west');    
        
        % write outputs
        geotiffwrite(fullfile(par.out_path,outdirs{ii},[outdirs{ii}, '.vel.geo.tif']),vel_regrid(y_ind,x_ind,ii),georef)
        geotiffwrite(fullfile(par.out_path,outdirs{ii},[outdirs{ii}, '.vstd.geo.tif']),vstd_regrid(y_ind,x_ind,ii),georef)
        geotiffwrite(fullfile(par.out_path,outdirs{ii},[outdirs{ii}, '.E.geo.tif']),compE_regrid(y_ind,x_ind,ii),georef)
        geotiffwrite(fullfile(par.out_path,outdirs{ii},[outdirs{ii}, '.N.geo.tif']),compN_regrid(y_ind,x_ind,ii),georef)
        geotiffwrite(fullfile(par.out_path,outdirs{ii},[outdirs{ii}, '.U.geo.tif']),compU_regrid(y_ind,x_ind,ii),georef)       
        
        disp([outdirs{ii} ' saved.'])
        
    end
    
end

if par.save_hgt == 1
    disp([par.out_path 'hgt.grd'])
    hgt(find(hgt==0)) = NaN;
    hgt = nanmedian(hgt,3);
    grdwrite2(x_regrid,y_regrid,hgt,[par.out_path 'hgt.grd']);
end
%%
plot_invLOS = 0;

if plot_invLOS == 1 && par.merge_tracks_along ~= 2 && par.ref2gnss ~= 0;
    for ii=1:nframes
        invlos = m_east .* compE_regrid(:,:,ii) + m_up .* compU_regrid(:,:,ii) + gnss_N .* compN_regrid(:,:,ii);
        invlos(find(mask_regrid(:,:,ii) == 0)) = NaN;
        figure()
        plt_data(x_regrid,y_regrid,vel_regrid(:,:,ii) - invlos,lonlim,latlim,[-5 5],'LOS residual (mm/yr)',fault_trace,borders)
        colormap(cpt.vik)
    end
end

%% end

disp('Run complete')
toc

end
