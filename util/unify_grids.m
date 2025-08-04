function [x_regrid,y_regrid,xx_regrid,yy_regrid,vel_regrid,vstd_regrid,...
    mask_regrid,compE_regrid,compU_regrid,compN_regrid,gnss_E,gnss_N,gnss_sE,gnss_sN, gnss_U, gnss_sU, hgt_regrid] ...
    = unify_grids(par,lon,lat,lon_comp,lat_comp,dx,dy,vel,vstd,mask,compE,compU,compN,gnss, hgt)
%=================================================================
% function unify_grids()
%-----------------------------------------------------------------
% Interpolate all inputs onto a common coordinate grid.
%                                                                  
% INPUT:                                                           
%   temp:
% OUTPUT:
%   temp:
%
% Andrew Watson     01-05-2023
%                                                                  
%=================================================================

nframes = length(vel);

%% define new grid

% limits and intervals for new grid
if par.auto_regrid == 1 % automate grid
    x_regrid = min(cellfun(@min,lon)) : min(cellfun(@min,dx)) : max(cellfun(@max,lon));
    y_regrid = min(cellfun(@min,lat)) : min(cellfun(@min,dy)) : max(cellfun(@max,lat));

elseif par.auto_regrid == 0   
    x_regrid = par.regrid_x;
    y_regrid = par.regrid_y;
    dx_regrid = x_regrid(2)-x_regrid(1); dy_regrid = y_regrid(2)-y_regrid(1);
    
    % check if requested grid spacing is near to input vels
    dx_check = min(cellfun(@min,dx)); dy_check = min(cellfun(@min,dy));
    if (abs(dx_check) - abs(dx_regrid)) > (dx_regrid*0.1) || ...
            (abs(dy_check) - abs(dy_regrid)) > (dy_regrid*0.1)
        warning('Resolution of input data is different from requested grid')
        disp(['Smallest input spacing in x = ' num2str(dx_check) ...
            ', requested spacing in x = ' num2str(dx_regrid)]);
        disp(['Smallest input spacing in y = ' num2str(dy_check) ...
            ', requested spacing in y = ' num2str(dy_regrid)]);
    end    
end

% new grid to project data onto
[xx_regrid,yy_regrid] = meshgrid(x_regrid,y_regrid); 

%% interpolate frames onto new grid

% pre-allocate
vel_regrid = single(zeros([size(xx_regrid) nframes]));
vstd_regrid = single(zeros(size(vel_regrid))); mask_regrid = single(zeros(size(vel_regrid)));
compE_regrid = single(zeros(size(vel_regrid))); compN_regrid = single(zeros(size(vel_regrid)));
compU_regrid = single(zeros(size(vel_regrid)));
hgt_regrid = single(zeros(size(vel_regrid)));

% for ii = 1:nframes
% 
%     % reduce interpolation area using coords of frame
%     [~,xind_min] = min(abs(x_regrid-lon{ii}(1)));
%     [~,xind_max] = min(abs(x_regrid-lon{ii}(end)));
%     [~,yind_min] = min(abs(y_regrid-lat{ii}(end)));
%     [~,yind_max] = min(abs(y_regrid-lat{ii}(1)));    
% 
%     % interpolate vel, vstd, and optional mask
%     [xx,yy] = meshgrid(lon{ii},lat{ii});
%     vel_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%         = interp2(xx,yy,vel{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%         yy_regrid(yind_min:yind_max,xind_min:xind_max));
% 
%     vstd_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%         = interp2(xx,yy,vstd{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%         yy_regrid(yind_min:yind_max,xind_min:xind_max));
% 
%     if par.use_mask == 1 || par.use_mask == 3
%         mask_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%             = interp2(xx,yy,mask{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%             yy_regrid(yind_min:yind_max,xind_min:xind_max));
%     end
% 
%     % ENU may not be downsampled by licsbas, hence different xx yy grids
%     [xx,yy] = meshgrid(lon_comp{ii},lat_comp{ii});
%     compE_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%         = interp2(xx,yy,compE{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%         yy_regrid(yind_min:yind_max,xind_min:xind_max));
%     compN_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%         = interp2(xx,yy,compN{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%         yy_regrid(yind_min:yind_max,xind_min:xind_max));
%     compU_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%         = interp2(xx,yy,compU{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%         yy_regrid(yind_min:yind_max,xind_min:xind_max));
%     if par.save_hgt == 1 || par.remove_linear_APS > 0
%         hgt_regrid(yind_min:yind_max,xind_min:xind_max,ii) ...
%             = interp2(xx,yy,hgt{ii},xx_regrid(yind_min:yind_max,xind_min:xind_max),...
%             yy_regrid(yind_min:yind_max,xind_min:xind_max));
%     end
% 
%     % report progress
%     if (mod(ii,round(nframes./10))) == 0
%         disp([num2str(round((ii./nframes)*100)) '% completed']);
%     end
% 
% end

for ii = 1:nframes

    % Get original lon/lat and data grids
    [xx, yy] = meshgrid(lon{ii}, lat{ii});
    x_flat = xx(:);
    y_flat = yy(:);
    
    % Flatten all data arrays
    vel_flat  = vel{ii}(:);
    vstd_flat = vstd{ii}(:);
    if par.use_mask == 1 || par.use_mask == 3
        mask_flat = mask{ii}(:);
    end

    % Find closest indices on regrid for each pixel
    [~, xi] = min(abs(x_regrid' - x_flat'), [], 1);  % size: [1 x N]
    [~, yi] = min(abs(y_regrid' - y_flat'), [], 1);  % size: [1 x N]

    % Loop over valid (non-NaN) values
    for k = 1:length(vel_flat)
        if ~isnan(vel_flat(k))
            vel_regrid(yi(k), xi(k), ii) = vel_flat(k);
            vstd_regrid(yi(k), xi(k), ii) = vstd_flat(k);
            if par.use_mask == 1 || par.use_mask == 3
                mask_regrid(yi(k), xi(k), ii) = mask_flat(k);
            end
        end
    end

    % ENU component: use lon_comp/lat_comp grid
    [xxc, yyc] = meshgrid(lon_comp{ii}, lat_comp{ii});
    xc_flat = xxc(:);
    yc_flat = yyc(:);
    compE_flat = compE{ii}(:);
    compN_flat = compN{ii}(:);
    compU_flat = compU{ii}(:);
    if par.save_hgt == 1 || par.remove_linear_APS > 0
        hgt_flat = hgt{ii}(:);
    end

    [~, xci] = min(abs(x_regrid' - xc_flat'), [], 1);
    [~, yci] = min(abs(y_regrid' - yc_flat'), [], 1);

    for k = 1:length(compE_flat)
        if ~isnan(compE_flat(k))
            compE_regrid(yci(k), xci(k), ii) = compE_flat(k);
            compN_regrid(yci(k), xci(k), ii) = compN_flat(k);
            compU_regrid(yci(k), xci(k), ii) = compU_flat(k);
            if par.save_hgt == 1 || par.remove_linear_APS > 0
                hgt_regrid(yci(k), xci(k), ii) = hgt_flat(k);
            end
        end
    end

    % Progress
    if mod(ii, round(nframes / 10)) == 0
        disp([num2str(round(ii / nframes * 100)), '% completed']);
    end

end


% crop regrid coords to area with valid pixels
if par.auto_regrid == 0 && par.crop_post_regrid == 1
    
    % new inds
    vel_crop = single(any(vel_regrid,3)); vel_crop(vel_crop==0) = nan;
    [~,x_ind,y_ind,~,~] = crop_nans(vel_crop,x_regrid,y_regrid);
    clear vel_crop
    
    % crop arrays
    vel_regrid = vel_regrid(y_ind,x_ind,:);
    vstd_regrid = vstd_regrid(y_ind,x_ind,:);
    mask_regrid = mask_regrid(y_ind,x_ind,:);
    compE_regrid = compE_regrid(y_ind,x_ind,:);
    compN_regrid = compN_regrid(y_ind,x_ind,:);
    compU_regrid = compU_regrid(y_ind,x_ind,:);
    x_regrid = x_regrid(x_ind);
    y_regrid = y_regrid(y_ind);
    xx_regrid = xx_regrid(y_ind,x_ind);
    yy_regrid = yy_regrid(y_ind,x_ind);
end

%% Interpolate GNSS
if isfield(gnss, 'x')
    if isvector(gnss.x) && isvector(gnss.y) && length(gnss.x) == length(gnss.y)
        % === Scattered GNSS points (1D vector) ===
        F_E = scatteredInterpolant(gnss.x, gnss.y, gnss.E, 'linear', 'none');
        F_N = scatteredInterpolant(gnss.x, gnss.y, gnss.N, 'linear', 'none');
        gnss_E = F_E(xx_regrid, yy_regrid);
        gnss_N = F_N(xx_regrid, yy_regrid);
        
        if isfield(gnss, 'U')
            F_U = scatteredInterpolant(gnss.x, gnss.y, gnss.U, 'linear', 'none');
            gnss_U = F_U(xx_regrid, yy_regrid);
        else
            gnss_U = [];
        end

        if par.gnss_uncer == 1
            F_sE = scatteredInterpolant(gnss.x, gnss.y, gnss.sE, 'linear', 1e6);
            F_sN = scatteredInterpolant(gnss.x, gnss.y, gnss.sN, 'linear', 1e6);
            gnss_sE = F_sE(xx_regrid, yy_regrid);
            gnss_sN = F_sN(xx_regrid, yy_regrid);

            if isfield(gnss, 'sU')
                F_sU = scatteredInterpolant(gnss.x, gnss.y, gnss.sU, 'linear', 1e6);
                gnss_sU = F_sU(xx_regrid, yy_regrid);
            else
                gnss_sU = [];
            end
        else
            gnss_sE = [];
            gnss_sN = [];
            gnss_sU = [];
        end

    else
        % === Gridded GNSS points (2D) ===
        [xx_gnss, yy_gnss] = meshgrid(gnss.x, gnss.y);
        gnss_E = interp2(xx_gnss, yy_gnss, gnss.E, xx_regrid, yy_regrid, 'linear', 0);
        gnss_N = interp2(xx_gnss, yy_gnss, gnss.N, xx_regrid, yy_regrid, 'linear', 0);

        if isfield(gnss, 'U')
            gnss_U = interp2(xx_gnss, yy_gnss, gnss.U, xx_regrid, yy_regrid, 'linear', 0);
        else
            gnss_U = [];
        end

        if par.gnss_uncer == 1
            gnss_sE = interp2(xx_gnss, yy_gnss, gnss.sE, xx_regrid, yy_regrid, 'linear', 1e6);
            gnss_sN = interp2(xx_gnss, yy_gnss, gnss.sN, xx_regrid, yy_regrid, 'linear', 1e6);
            if isfield(gnss, 'sU')
                gnss_sU = interp2(xx_gnss, yy_gnss, gnss.sU, xx_regrid, yy_regrid, 'linear', 1e6);
            else
                gnss_sU = [];
            end
        else
            gnss_sE = [];
            gnss_sN = [];
            gnss_sU = [];
        end
    end
else
    % Empty GNSS case
    gnss_E = [];
    gnss_N = [];
    gnss_U = [];
    gnss_sE = [];
    gnss_sN = [];
    gnss_sU = [];
end

%% format

% change zeros to nans
vstd_regrid(vstd_regrid==0) = nan;
compE_regrid(compE_regrid==0) = nan;
compN_regrid(compN_regrid==0) = nan;
compU_regrid(compU_regrid==0) = nan;


end